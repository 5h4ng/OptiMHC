import logging
from collections import defaultdict
from typing import Dict, List, Tuple, Union

import networkx as nx
import numpy as np
import pandas as pd
from scipy.stats import entropy
from tqdm import tqdm

from optimhc import utils
from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.feature.factory import feature_generator_factory
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


class OverlappingPeptideFeatureGenerator(BaseFeatureGenerator):
    """
    Generates features based on peptide sequence overlaps using the Overlap-Layout-Consensus (OLC) algorithm.

    This generator constructs an overlap graph of peptides, removes transitive edges, simplifies the graph to contigs,
    and computes features such as the number of overlaps, log-transformed overlap counts, overlap ranks, and contig lengths.
    It also filters out peptides with low entropy or outlier lengths before processing.
    Additionally, it records detailed information about brother peptides and contigs, accessible via the get_all_data method.

    Parameters
    ----------
    peptides : list of str
        List of peptide sequences.
    min_overlap_length : int, optional
        Minimum required overlap length for peptides to be considered overlapping. Default is 6.
    min_length : int, optional
        Minimum peptide length to include in processing. Default is 7.
    max_length : int, optional
        Maximum peptide length to include in processing. Default is 60.
    min_entropy : float, optional
        Minimum Shannon entropy for peptides to include in processing. Default is 0.
    fill_missing : str, optional
        Method to fill missing values for filtered peptides. Options are 'median' or 'zero'. Default is 'median'.
    remove_pre_nxt_aa : bool, optional
        Whether to remove the preceding and following amino acids from peptides. Default is False.
    remove_modification : bool, optional
        Whether to remove modifications from peptides. Default is True.

    Attributes
    ----------
    original_peptides : list of str
        Original list of peptide sequences.
    min_overlap_length : int
        Minimum required overlap length.
    min_length : int
        Minimum peptide length.
    max_length : int
        Maximum peptide length.
    min_entropy : float
        Minimum Shannon entropy.
    fill_missing : str
        Method to fill missing values.
    remove_pre_nxt_aa : bool
        Whether to remove preceding and following amino acids.
    remove_modification : bool
        Whether to remove modifications.
    filtered_peptides : list of str
        List of peptides after filtering.
    filtered_indices : list of int
        Indices of filtered peptides.
    peptide_to_index : dict of str to int
        Mapping of peptides to their indices.
    overlap_data : pd.DataFrame
        DataFrame containing overlap data.
    peptide_to_contig : dict of str to int
        Mapping of peptides to their contig indices.
    assembled_contigs : list of dict
        List of assembled contigs.
    full_data : pd.DataFrame
        Full data including brother peptides and contig information.
    _overlap_graph : nx.DiGraph
        Overlap graph.
    _simplified_graph : nx.DiGraph
        Simplified graph with transitive edges removed.

    Notes
    -----
    Key Data Structures:
        1. contigs: List[List[str]]
           - Represents non-branching paths in the overlap graph
           - Each inner list contains peptide sequences that form a continuous chain
           - Example: [['PEPTIDE1', 'PEPTIDE2'], ['PEPTIDE3']]

        2. assembled_contigs: List[Dict]
           - Contains the assembled sequences and their constituent peptides
           - Each dictionary has two keys:
               'sequence': The merged/assembled sequence of overlapping peptides
               'peptides': List of peptides that were used to build this contig
           - Example: [
               {
                   'sequence': 'LONGPEPTIDESEQUENCE',
                   'peptides': ['LONGPEP', 'PEPTIDE', 'SEQUENCE']
               },
               {
                   'sequence': 'SINGLEPEPTIDE',
                   'peptides': ['SINGLEPEPTIDE']
               }
           ]

        3. peptide_to_contig: Dict[str, int]
           - Maps each peptide to its contig index in assembled_contigs
           - Key: peptide sequence
           - Value: index of the contig containing this peptide
           - Example: {
               'LONGPEP': 0,
               'PEPTIDE': 0,
               'SEQUENCE': 0,
               'SINGLEPEPTIDE': 1
           }

        4. overlap_graph (_overlap_graph): nx.DiGraph
           - Directed graph representing all possible overlaps between peptides
           - Nodes: peptide sequences
           - Edges: overlaps between peptides
           - Edge weights: length of overlap

        5. simplified_graph (_simplified_graph): nx.DiGraph
           - Simplified version of overlap_graph with transitive edges removed
           - Used for final contig assembly
           - More efficient representation of essential overlaps
    """

    def __init__(
        self,
        peptides: List[str],
        min_overlap_length: int = 6,
        min_length: int = 7,
        max_length: int = 60,
        min_entropy: float = 0,
        fill_missing: str = "median",  # 'median' or 'zero'
        remove_pre_nxt_aa: bool = False,
        remove_modification: bool = True,
        *args,
        **kwargs,
    ):
        self.original_peptides = peptides
        self.min_overlap_length = min_overlap_length
        self.min_length = min_length
        self.max_length = max_length
        self.min_entropy = min_entropy
        self.fill_missing = fill_missing.lower()
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.remove_modification = remove_modification
        self.filtered_peptides = []
        self.filtered_indices = []
        self.peptide_to_index = {}
        self.overlap_data = None
        self.peptide_to_contig = {}
        self.assembled_contigs = []
        self.full_data = None
        self._overlap_graph = None
        self._simplified_graph = None
        logger.info(
            f"Initialized OverlappingPeptideFeatureGenerator with {len(peptides)} peptides and minimum overlap length: {min_overlap_length}"
        )
        logger.info(
            f"remove_pre_nxt_aa: {remove_pre_nxt_aa}, remove_modification: {remove_modification}"
        )
        logger.info(
            f"Peptide filtering parameters - min_length: {min_length}, max_length: {max_length}, min_entropy: {min_entropy}"
        )

    @property
    def id_column(self) -> List[str]:
        """
        Returns a list of input columns required for the feature generator.

        Returns:
            List[str]: List of input columns.
        """
        return ["Peptide"]

    @property
    def feature_columns(self) -> List[str]:
        """Returns the feature column names."""
        """
        return ['contig_member_count', 'log_contig_member_count', 'contig_member_rank', 'log_contig_member_rank',
                 'contig_seq_length_diff', 'contig_length', 'contig_extension_ratio']
        """
        return [
            "contig_member_count",
            "contig_extension_ratio",
            "contig_member_rank",
            "contig_length",
        ]

    @property
    def overlap_graph(self) -> nx.DiGraph:
        """Returns the overlap graph."""
        return self._overlap_graph

    @property
    def simplified_graph(self) -> nx.DiGraph:
        """Returns the layout graph."""
        return self._simplified_graph

    @property
    def contigs(self) -> List[Dict]:
        """Returns the assembled contigs."""
        return self.assembled_contigs

    def _shannon_entropy(self, sequence: str) -> float:
        """
        Calculate the Shannon entropy of a peptide sequence.

        Parameters:
            sequence (str): Peptide sequence.

        Returns:
            float: Shannon entropy value.
        """
        bases = list(set(sequence))
        freq_list = [sequence.count(base) / len(sequence) for base in bases]
        return entropy(freq_list, base=2)

    def _preprocess_peptides(self, peptide: str) -> str:
        if self.remove_pre_nxt_aa:
            peptide = utils.strip_flanking_and_charge(peptide)
        if self.remove_modification:
            peptide = utils.remove_modifications(peptide)
        # U -> C
        peptide = peptide.replace("U", "C")
        return peptide

    def _filter_peptides(self, peptides: List[str]) -> List[str]:
        """
        Filter out peptides based on length and entropy.

        Parameters:
            peptides (List[str]): List of peptide sequences.

        Returns:
            List[str]: Filtered list of peptide sequences.
        """
        filtered_peptides = []
        for peptide in peptides:
            if len(peptide) < self.min_length or len(peptide) > self.max_length:
                continue
            entropy_val = self._shannon_entropy(peptide)
            if entropy_val < self.min_entropy:
                continue
            filtered_peptides.append(peptide)
        logger.info(
            f"Filtered out {len(peptides) - len(filtered_peptides)} peptides based on length and entropy."
        )
        logger.info(f"Remaining peptides: {len(filtered_peptides)}")
        return filtered_peptides

    def _construct_prefix_index(
        self, peptides: List[str], min_overlap_length: int
    ) -> Dict[str, List[int]]:
        """
        Construct an index of prefixes for all peptides.

        Parameters:
            peptides (List[str]): List of peptide sequences.

        Returns:
            Dict[str, List[int]]: Dictionary mapping prefixes to list of peptide indices.
        """
        prefix_index = defaultdict(list)
        for idx, seq in enumerate(peptides):
            seq_len = len(seq)
            # Store prefixes of length from min_overlap_length up to seq_len
            for i in range(min_overlap_length, seq_len + 1):
                prefix = seq[:i]
                prefix_index[prefix].append(idx)
        return prefix_index

    def _build_overlap_graph(
        self, peptides: List[str], prefix_index: Dict[str, List[int]]
    ) -> nx.DiGraph:
        """
        Build the overlap graph from the list of peptides.

        Parameters:
            peptides (List[str]): List of peptide sequences.
            prefix_index (Dict[str, List[int]]): Index of prefixes.

        Returns:
            nx.DiGraph: Overlap graph.
        """
        G = nx.DiGraph()
        for idx, peptide in enumerate(peptides):
            seq_len = len(peptide)
            G.add_node(peptide)
            # Iterate over possible suffix lengths
            for i in range(self.min_overlap_length, seq_len):
                suffix = peptide[-i:]
                if suffix in prefix_index:
                    for matching_idx in prefix_index[suffix]:
                        if matching_idx != idx:
                            matching_peptide = peptides[matching_idx]
                            overlap_length = len(suffix)
                            # Ensure the overlap is the maximal possible
                            existing_weight = (
                                G[peptide][matching_peptide]["weight"]
                                if G.has_edge(peptide, matching_peptide)
                                else 0
                            )
                            if overlap_length > existing_weight:
                                G.add_edge(peptide, matching_peptide, weight=overlap_length)
            # Handle the case where the entire peptide matches the prefix of another peptide
            suffix = peptide  # Full peptide as suffix
            if suffix in prefix_index:
                for matching_idx in prefix_index[suffix]:
                    if matching_idx != idx:
                        matching_peptide = peptides[matching_idx]
                        G.add_edge(peptide, matching_peptide, weight=seq_len)
        return G

    def _remove_transitive_edges(self, G: nx.DiGraph) -> nx.DiGraph:
        """
        Remove transitive edges from the overlap graph G.

        Parameters:
            G (nx.DiGraph): Overlap graph.

        Returns:
            nx.DiGraph: Simplified graph with transitive edges removed.
        """
        logger.info("Removing transitive edges from the overlap graph.")

        ## if A->B->C, A->C, then remove A->C
        to_remove = []
        for u in G:
            for v in G.successors(u):
                for w in G.successors(v):
                    if G.has_edge(u, w):
                        to_remove.append((u, w))
        G.remove_edges_from(to_remove)

        to_remove = []
        for u in G:
            for v in G.successors(u):
                for w in G.successors(v):
                    for x in G.successors(w):
                        if G.has_edge(u, x):
                            to_remove.append((u, x))
        G.remove_edges_from(to_remove)

        return G

    def _simplify_graph_to_contigs(self, G: nx.DiGraph) -> List[List[str]]:
        """
        Simplify the graph by finding non-branching paths (contigs).

        Parameters:
            G (nx.DiGraph): Overlap graph.

        Returns:
            List[List[str]]: List of contigs, each contig is a list of peptide sequences.
        """
        logger.info("Simplifying graph to contigs (non-branching paths).")
        contigs = []
        visited = set()

        for node in G.nodes():
            # Check if node is a potential starting point of a contig
            if G.in_degree(node) != 1 or G.out_degree(node) != 1:
                if node not in visited:
                    contig = [node]
                    visited.add(node)
                    # Extend contig forward
                    current_node = node
                    while G.out_degree(current_node) == 1:
                        successor = next(G.successors(current_node))
                        if G.in_degree(successor) == 1 and successor not in visited:
                            contig.append(successor)
                            visited.add(successor)
                            current_node = successor
                        else:
                            break
                    contigs.append(contig)

        logger.info(f"Found {len(contigs)} contigs.")
        return contigs

    def _assemble_contigs(self, contigs: List[List[str]], G: nx.DiGraph) -> List[Dict]:
        """
        Assemble sequences for each contig.

        Parameters:
            contigs (List[List[str]]): List of contigs.

        Returns:
            List[Dict]: List of assembled contigs with sequence and peptides.
        """
        logger.info("Assembling contigs from peptides.")
        assembled_contigs = []
        for contig in contigs:
            if len(contig) == 1:
                # Single node contig
                assembled_seq = contig[0]
            else:
                assembled_seq = contig[0]
                for i in range(1, len(contig)):
                    # Get the overlap length between consecutive peptides
                    overlap_length = G[contig[i - 1]][contig[i]]["weight"]
                    # Append non-overlapping part of the next peptide
                    assembled_seq += contig[i][overlap_length:]
            assembled_contigs.append({"sequence": assembled_seq, "peptides": contig})
        return assembled_contigs

    def _map_peptides_to_contigs(self, assembled_contigs: List[Dict]):
        """
        Map each peptide to its contig.

        Parameters:
            assembled_contigs (List[Dict]): List of assembled contigs.
        """
        logger.info("Mapping peptides to contigs.")
        peptide_to_contig = {}
        for idx, contig in enumerate(assembled_contigs):
            peptides_in_contig = contig["peptides"]
            for peptide in peptides_in_contig:
                peptide_to_contig[peptide] = idx
        return peptide_to_contig

    def _remove_redundant_peptides(self, peptides: List[str]) -> Tuple[List[str], Dict[str, str]]:
        """
        Remove peptides that are fully contained in other peptides.

        Parameters:
            peptides (List[str]): List of peptide sequences.

        Returns:
            accepted_peptides (List[str]): List of accepted peptides not redundant (largest container peptides).
            redundant_mapping (Dict[str, str]): Mapping of redundant peptide to its largest container peptide.
        """
        sorted_peptides = sorted(peptides, key=len, reverse=True)
        min_pep_len = len(sorted_peptides[-1])
        peptides_set = set(sorted_peptides)
        accepted_set = set()
        redundant_mapping = {}
        for pep in sorted_peptides:
            if pep not in peptides_set:
                continue
            if pep in accepted_set:
                continue
            accepted_set.add(pep)
            pep_len = len(pep)
            for L in range(min_pep_len, pep_len + 1):
                for i in range(0, pep_len - L + 1):
                    sub_pep = pep[i : i + L]
                    if sub_pep == pep:
                        continue
                    if sub_pep in peptides_set:
                        redundant_mapping[sub_pep] = pep
                        peptides_set.remove(sub_pep)

        logger.info(f"Remove {len(redundant_mapping)} redundant peptides.")
        logger.info(f"Accepted {len(accepted_set)} non-redundant peptides.")
        return list(accepted_set), redundant_mapping

    def _build_full_contig_map(self, peptides: List[str]) -> Dict[int, List[str]]:
        """
        Build a mapping from contig index to list of peptides (both accepted and redundant).

        Parameters:
            peptides (List[str]): Original list of peptides.

        Returns:
            Dict[int, List[str]]: Mapping from contig index to list of peptides.
        """
        full_contig_map = defaultdict(list)
        for pep in peptides:
            contig_idx = self.peptide_to_contig.get(pep, None)
            if contig_idx is not None:
                full_contig_map[contig_idx].append(pep)
        return full_contig_map

    def _calculate_overlap_contig_features(self, peptides: List[str]) -> pd.DataFrame:
        """
        Calculate overlap and contig related features for peptides.

        This method computes the overlap graph, simplifies it to contigs, assembles contigs,
        and generates features for each peptide. It also handles redundant peptides by mapping them
        to their container peptides.

        Parameters:
            peptides (List[str]): List of peptide sequences.

        Returns:
            pd.DataFrame: DataFrame containing the computed features.
        """
        accepted_peptides, redundant_mapping = self._remove_redundant_peptides(peptides)

        logger.info("Constructing prefix index...")
        prefix_index = self._construct_prefix_index(accepted_peptides, self.min_overlap_length)

        logger.info("Building overlap graph...")
        self._overlap_graph = self._build_overlap_graph(accepted_peptides, prefix_index)

        logger.info(
            f"Overlap graph has {self._overlap_graph.number_of_nodes()} nodes and {self._overlap_graph.number_of_edges()} edges."
        )
        self._simplified_graph = self._remove_transitive_edges(self._overlap_graph)

        logger.info(
            f"Simplified graph has {self._simplified_graph.number_of_nodes()} nodes and {self._simplified_graph.number_of_edges()} edges."
        )
        contigs = self._simplify_graph_to_contigs(self._simplified_graph)

        logger.info(f"Found {len(contigs)} contigs.")
        self.assembled_contigs = self._assemble_contigs(contigs, self._simplified_graph)
        self.peptide_to_contig = self._map_peptides_to_contigs(self.assembled_contigs)

        # Map redundant peptides to their container peptides
        for redundant_peptide, container_peptide in redundant_mapping.items():
            if container_peptide not in self.peptide_to_contig:
                logger.debug(f"Container peptide {container_peptide} not found in contigs.")
                logger.debug(
                    "This may occur if the container peptide is a branching node in the overlap graph."
                )
                logger.debug(f"Assigning {container_peptide} to its own contig.")

                new_contig_idx = len(self.assembled_contigs)
                new_contig = {
                    "sequence": container_peptide,
                    "peptides": [container_peptide],
                    "full_contig_peptides": [container_peptide],
                }
                self.assembled_contigs.append(new_contig)
                self.peptide_to_contig[container_peptide] = new_contig_idx

            self.peptide_to_contig[redundant_peptide] = self.peptide_to_contig[container_peptide]

        # Build a mapping from contig index to list of peptides (both accepted and redundant)
        full_contig_map = self._build_full_contig_map(peptides)

        for contig_index, full_peptides in full_contig_map.items():
            self.assembled_contigs[contig_index]["full_contig_peptides"] = full_peptides

        # Compute features for each peptide in the filtered list
        feature_list = []
        for pep in peptides:
            contig_idx = self.peptide_to_contig.get(pep, None)
            if contig_idx is not None:
                full_count = len(full_contig_map[contig_idx])
                contig_member_count = full_count
                contig_length = len(self.assembled_contigs[contig_idx]["sequence"])
            else:
                contig_member_count = 0
                contig_length = len(pep)
            feature_list.append(
                {
                    "clean_peptide": pep,
                    "contig_member_count": contig_member_count,
                    "contig_length": contig_length,
                }
            )

        features_df = pd.DataFrame(feature_list)
        features_df["log_contig_member_count"] = features_df["contig_member_count"].apply(
            lambda x: np.log(x + 1e-6)
        )
        features_df["contig_member_rank"] = features_df["contig_member_count"].rank(
            method="min", ascending=False
        )
        features_df["log_contig_member_rank"] = features_df["contig_member_rank"].apply(
            lambda x: np.log(x + 1e-6)
        )
        features_df["contig_seq_length_diff"] = features_df["contig_length"] - features_df[
            "clean_peptide"
        ].apply(len)
        features_df["contig_extension_ratio"] = features_df[
            "contig_seq_length_diff"
        ] / features_df["clean_peptide"].apply(len)

        return features_df

    def _integrate_overlap_features(self) -> pd.DataFrame:
        """
        Compute the features for each peptide, ensuring output aligns with input.

        This method preprocesses and filters peptides, calculates overlap and contig features,
        maps them back to the original peptides, and fills missing values.

        Returns:
            pd.DataFrame: DataFrame containing features for each peptide.
        """
        if self.overlap_data is None:
            # 1. Preprocess and filter peptides
            self.overlap_data = pd.DataFrame(self.original_peptides, columns=["Peptide"])
            self.overlap_data["clean_peptide"] = self.overlap_data["Peptide"].apply(
                self._preprocess_peptides
            )
            self.filtered_peptides = self._filter_peptides(
                self.overlap_data["clean_peptide"].unique().tolist()
            )

            # 2. Compute overlaps and features for filtered peptides
            features_df = self._calculate_overlap_contig_features(self.filtered_peptides)

            # 3. Map features back to the original peptides
            logger.info("Mapping features back to original peptides.")
            self.overlap_data = self.overlap_data.merge(
                features_df, on="clean_peptide", how="left"
            )

            # 4. Handle missing features for peptides that were filtered out
            missing_counts = self.overlap_data["contig_member_count"].isna().sum()
            logger.info(
                f"Number of peptides with missing features (filtered out): {missing_counts}"
            )
            if self.fill_missing == "median":
                logger.info("Filling missing values with median.")
                median_values = {
                    "contig_member_count": self.overlap_data["contig_member_count"].median(),
                    "log_contig_member_count": self.overlap_data[
                        "log_contig_member_count"
                    ].median(),
                    "contig_member_rank": self.overlap_data["contig_member_rank"].median(),
                    "log_contig_member_rank": self.overlap_data["log_contig_member_rank"].median(),
                    "contig_length": self.overlap_data["contig_length"].median(),
                    "contig_seq_length_diff": self.overlap_data["contig_seq_length_diff"].median(),
                    "contig_extension_ratio": self.overlap_data["contig_extension_ratio"].median(),
                }
                self.overlap_data.fillna(value=median_values, inplace=True)
            elif self.fill_missing == "zero":
                logger.info("Filling missing values with zero.")
                self.overlap_data.fillna(value=0, inplace=True)
            else:
                logger.warning(
                    f"Invalid fill_missing option '{self.fill_missing}'. Defaulting to zero."
                )
                self.overlap_data.fillna(value=0, inplace=True)
            logger.info("Feature computation completed.")
        else:
            logger.info("Features have already been computed. Skipping recomputation.")
        return self.overlap_data

    def generate_features(self) -> pd.DataFrame:
        """
        Generates features for peptide overlaps, including the count of overlapping peptides, contig length,
        and log-transformed counts and ranks.

        Returns:
            pd.DataFrame: DataFrame containing the features.
        """
        features_df = self._integrate_overlap_features()
        features_df = features_df[["Peptide"] + self.feature_columns]
        logger.info(f"Generated overlap features for {len(features_df)} peptides.")
        return features_df

    def get_full_data(self) -> pd.DataFrame:
        """
        Returns the full data including brother peptides and contig information for each peptide.
        In the output, the lists of contig peptides and brother peptides include redundant peptides,
        so that their counts match the corresponding peptide and contig_member_count.

        Returns:
            pd.DataFrame: DataFrame containing peptides and their brother peptides and contigs.
        """
        self._integrate_overlap_features()
        if self.full_data is not None:
            logger.info("Full data has already been computed. Returning cached data.")
            return self.full_data
        data_list = []

        for peptide in tqdm(self.filtered_peptides):
            contig_idx = self.peptide_to_contig.get(peptide, None)
            if contig_idx is not None:
                contig_info = self.assembled_contigs[contig_idx]
                # Use full contig peptides (including redundant ones) if available
                full_peptides = contig_info.get("full_contig_peptides", contig_info["peptides"])
                brother_peptides = [p for p in full_peptides if p != peptide]
                data_list.append(
                    {
                        "clean_peptide": peptide,
                        "BrotherPeptides": brother_peptides,
                        "ContigSequence": contig_info["sequence"],
                        "ContigPeptides": full_peptides,
                    }
                )

        full_data_df = pd.DataFrame(data_list)
        self.full_data = self.overlap_data.merge(full_data_df, on="clean_peptide", how="left")
        return self.full_data

    @classmethod
    def from_config(cls, psms, config, params):
        instance = cls(
            peptides=list(set(psms.peptides)),
            min_overlap_length=params.get("minOverlapLength", 8),
            min_length=params.get("minLength", 8),
            max_length=params.get("maxLength", 25),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
        )
        instance._overlapping_score = params.get("overlappingScore", None)
        return instance

    def apply(self, psms, source):
        features = self.generate_features()
        full_data = self.get_full_data()

        psms.add_metadata(
            full_data[["Peptide", "contig_member_count", "ContigSequence"]],
            psms_key=psms.peptide_column,
            metadata_key="Peptide",
            source=source,
        )
        psms.add_features(
            features,
            psms_key=psms.peptide_column,
            feature_key=self.id_column,
            source=source,
        )

        if self._overlapping_score:
            assign_brother_aggregated_feature(
                psms,
                feature_columns=self._overlapping_score,
                overlapping_source=source,
                source_name="ContigFeatures",
            )


def assign_brother_aggregated_feature(
    psms: PsmContainer,
    feature_columns: Union[str, List[str]],
    overlapping_source: str,
    source_name: str = "OverlappingGroupFeatures",
) -> None:
    """
    Assign aggregated features based on brother peptides to the PSMs.

    For PSMs with the same ContigSequence (brother peptides), compute the mean of specified features
    and assign these aggregated features back to each PSM in the group. Additionally, compute
    the sum as mean * (contig_member_count + 1). If a PSM does not have a ContigSequence (no brothers),
    its new features will be set to the original values.

    Parameters:
        psms (PsmContainer): PSM container containing the peptides and features.
        feature_columns (Union[str, List[str]]): Name of the feature column(s) to aggregate.
        overlapping_source (str): Source name of the overlapping peptide features.
        source_name (str): Name of the new feature source.

    Returns:
        None
    """
    if isinstance(feature_columns, str):
        feature_columns = [feature_columns]
    psms_df = psms.psms

    if psms.metadata_column is None:
        raise ValueError("The PSMs do not contain metadata.")
    metadata = psms_df[psms.metadata_column]

    def get_overlapping_data(x):
        if isinstance(x, dict):
            return x.get(overlapping_source, {})
        else:
            logger.warning(f"Invalid metadata entry: {x}")
            return {}

    overlapping_data = metadata.apply(get_overlapping_data)

    def get_contig_sequence(x):
        if isinstance(x, dict):
            return x.get("ContigSequence", None)
        else:
            logger.warning(f"Invalid overlapping data entry: {x}")
            return None

    contig_sequences = overlapping_data.apply(get_contig_sequence)

    psms_df["ContigSequence"] = contig_sequences

    if "contig_member_count" not in psms_df.columns:
        raise ValueError("'contig_member_count' column not found in PSMs.")

    missing_features = [feature for feature in feature_columns if feature not in psms_df.columns]
    if missing_features:
        raise ValueError(f"Feature columns not found in PSMs: {missing_features}")

    grouped_mean = psms_df.groupby("ContigSequence")[feature_columns].mean().reset_index()
    grouped_mean = grouped_mean.rename(
        columns={feature: f"{feature}_contig_avg" for feature in feature_columns}
    )

    psms_with_agg = psms_df.merge(grouped_mean, on="ContigSequence", how="left")

    for feature in feature_columns:
        mean_feature = f"{feature}_contig_avg"
        sum_feature = f"{feature}_contig_sum"
        psms_with_agg["contig_member_count"] = psms_with_agg["contig_member_count"].fillna(0)
        psms_with_agg[sum_feature] = (
            psms_with_agg[mean_feature] * (psms_with_agg["contig_member_count"])
        )
        psms_with_agg[sum_feature].fillna(psms_with_agg[feature], inplace=True)

    for feature in feature_columns:
        mean_feature = f"{feature}_contig_avg"
        psms_with_agg[mean_feature].fillna(psms_with_agg[feature], inplace=True)

    agg_feature_columns = [f"{feature}_contig_avg" for feature in feature_columns] + [
        f"{feature}_contig_sum" for feature in feature_columns
    ]

    new_features_df = psms_with_agg[agg_feature_columns]

    psms.add_features_by_index(features_df=new_features_df, source=source_name)


feature_generator_factory.register_generator(
    "OverlappingPeptide", OverlappingPeptideFeatureGenerator
)
