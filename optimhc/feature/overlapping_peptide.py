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
    """Generate layout features from exact suffix-prefix peptide overlaps.

    The generator filters and normalizes peptide sequences, removes contained
    sequences from layout construction, builds a directed overlap graph, performs
    offset-consistent transitive reduction, and spells deterministic non-branching
    layouts. Contained peptides are mapped back before features are calculated.

    This is an OLC-inspired overlap/layout heuristic. It does not perform a
    residue-level consensus step; assembled sequences are layout sequences.
    """

    def __init__(
        self,
        peptides: List[str],
        min_overlap_length: int = 6,
        min_length: int = 7,
        max_length: int = 60,
        min_entropy: float = 0,
        fill_missing: str = "median",
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
        """Return the input identifier column."""
        return ["Peptide"]

    @property
    def feature_columns(self) -> List[str]:
        """Return model feature columns."""
        return [
            "contig_member_count",
            "contig_extension_ratio",
            "contig_member_rank",
            "contig_length",
        ]

    @property
    def overlap_graph(self) -> nx.DiGraph:
        """Return the full overlap graph."""
        return self._overlap_graph

    @property
    def simplified_graph(self) -> nx.DiGraph:
        """Return the transitively reduced layout graph."""
        return self._simplified_graph

    @property
    def contigs(self) -> List[Dict]:
        """Return assembled layouts."""
        return self.assembled_contigs

    def _shannon_entropy(self, sequence: str) -> float:
        """Calculate sequence Shannon entropy in bits."""
        bases = list(set(sequence))
        freq_list = [sequence.count(base) / len(sequence) for base in bases]
        return entropy(freq_list, base=2)

    def _preprocess_peptides(self, peptide: str) -> str:
        if self.remove_pre_nxt_aa:
            peptide = utils.strip_flanking_and_charge(peptide)
        if self.remove_modification:
            peptide = utils.remove_modifications(peptide)
        # Preserve the legacy U/C equivalence used by overlap features.
        peptide = peptide.replace("U", "C")
        return peptide

    def _filter_peptides(self, peptides: List[str]) -> List[str]:
        """Return unique peptides passing the length and entropy filters."""
        filtered_peptides = []
        for peptide in sorted(set(peptides)):
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
        """Map eligible prefixes to peptide indexes."""
        prefix_index = defaultdict(list)
        for idx, seq in enumerate(peptides):
            seq_len = len(seq)
            for i in range(min_overlap_length, seq_len + 1):
                prefix = seq[:i]
                prefix_index[prefix].append(idx)
        return prefix_index

    def _build_overlap_graph(
        self, peptides: List[str], prefix_index: Dict[str, List[int]]
    ) -> nx.DiGraph:
        """Build the directed maximal suffix-prefix overlap graph."""
        G = nx.DiGraph()
        for idx, peptide in enumerate(peptides):
            seq_len = len(peptide)
            G.add_node(peptide)
            for i in range(self.min_overlap_length, seq_len):
                suffix = peptide[-i:]
                if suffix in prefix_index:
                    for matching_idx in prefix_index[suffix]:
                        if matching_idx != idx:
                            matching_peptide = peptides[matching_idx]
                            overlap_length = len(suffix)
                            existing_weight = (
                                G[peptide][matching_peptide]["weight"]
                                if G.has_edge(peptide, matching_peptide)
                                else 0
                            )
                            if overlap_length > existing_weight:
                                G.add_edge(peptide, matching_peptide, weight=overlap_length)
            suffix = peptide
            if suffix in prefix_index:
                for matching_idx in prefix_index[suffix]:
                    if matching_idx != idx:
                        matching_peptide = peptides[matching_idx]
                        G.add_edge(peptide, matching_peptide, weight=seq_len)
        return G

    def _remove_transitive_edges(self, G: nx.DiGraph) -> nx.DiGraph:
        """Remove edges implied by an alternative path at the same sequence offset."""
        logger.info("Removing transitive edges from the overlap graph.")
        reduced = G.copy()

        def edge_offset(source: str, edge_data: Dict) -> int:
            return len(source) - edge_data["weight"]

        def has_offset_consistent_path(source: str, target: str, target_offset: int) -> bool:
            stack = [(source, 0)]
            visited_states = {(source, 0)}

            while stack:
                node, current_offset = stack.pop()
                for successor in sorted(reduced.successors(node), reverse=True):
                    step = edge_offset(node, reduced[node][successor])
                    next_offset = current_offset + step
                    if next_offset > target_offset:
                        continue
                    if successor == target and next_offset == target_offset:
                        return True
                    state = (successor, next_offset)
                    if next_offset < target_offset and state not in visited_states:
                        visited_states.add(state)
                        stack.append(state)
            return False

        edges = sorted(
            reduced.edges(data=True),
            key=lambda edge: (edge_offset(edge[0], edge[2]), edge[0], edge[1]),
        )
        for source, target, edge_data in edges:
            if not reduced.has_edge(source, target):
                continue
            target_offset = edge_offset(source, edge_data)
            saved_data = dict(reduced[source][target])
            reduced.remove_edge(source, target)
            if not has_offset_consistent_path(source, target, target_offset):
                reduced.add_edge(source, target, **saved_data)

        return reduced

    def _simplify_graph_to_contigs(self, G: nx.DiGraph) -> List[List[str]]:
        """Extract deterministic non-branching layouts and remaining cycles."""
        logger.info("Simplifying graph to contigs (non-branching paths).")
        contigs = []
        visited = set()

        def has_unambiguous_incoming_edge(node: str) -> bool:
            if G.in_degree(node) != 1:
                return False
            predecessor = next(iter(G.predecessors(node)))
            return G.out_degree(predecessor) == 1

        # TODO: Consider branch-related problems
        # the grouping rule and features are still TBD.
        for node in sorted(G.nodes()):
            # An unambiguous incoming edge joins this node to its predecessor's
            # unitig, so this node must not start another one.
            if node in visited or has_unambiguous_incoming_edge(node):
                continue
            contig = [node]
            visited.add(node)
            current_node = node
            while G.out_degree(current_node) == 1:
                successor = next(iter(sorted(G.successors(current_node))))
                if G.in_degree(successor) == 1 and successor not in visited:
                    contig.append(successor)
                    visited.add(successor)
                    current_node = successor
                else:
                    break
            contigs.append(contig)

        # Components consisting only of 1-in/1-out nodes have no natural path head.
        # Emit each remaining cycle once, starting at its lexicographically smallest node.
        for node in sorted(G.nodes()):
            if node in visited:
                continue
            contig = []
            current_node = node
            while current_node not in visited:
                contig.append(current_node)
                visited.add(current_node)
                if G.out_degree(current_node) != 1:
                    break
                current_node = next(iter(sorted(G.successors(current_node))))
            if contig:
                contigs.append(contig)

        logger.info(f"Found {len(contigs)} contigs.")
        return contigs

    def _assemble_contigs(self, contigs: List[List[str]], G: nx.DiGraph) -> List[Dict]:
        """Spell each layout by appending non-overlapping peptide suffixes."""
        logger.info("Assembling contigs from peptides.")
        assembled_contigs = []
        for contig in contigs:
            if len(contig) == 1:
                assembled_seq = contig[0]
            else:
                assembled_seq = contig[0]
                for i in range(1, len(contig)):
                    overlap_length = G[contig[i - 1]][contig[i]]["weight"]
                    assembled_seq += contig[i][overlap_length:]
            assembled_contigs.append({"sequence": assembled_seq, "peptides": contig})
        return assembled_contigs

    def _map_peptides_to_contigs(self, assembled_contigs: List[Dict]):
        """Map each peptide to its assembled layout index."""
        logger.info("Mapping peptides to contigs.")
        peptide_to_contig = {}
        for idx, contig in enumerate(assembled_contigs):
            peptides_in_contig = contig["peptides"]
            for peptide in peptides_in_contig:
                peptide_to_contig[peptide] = idx
        return peptide_to_contig

    def _remove_redundant_peptides(self, peptides: List[str]) -> Tuple[List[str], Dict[str, str]]:
        """Remove contained sequences from layout and return their container mapping."""
        sorted_peptides = sorted(set(peptides), key=lambda peptide: (-len(peptide), peptide))
        if not sorted_peptides:
            return [], {}
        min_pep_len = len(sorted_peptides[-1])
        peptides_set = set(sorted_peptides)
        accepted_peptides = []
        redundant_mapping = {}
        for pep in sorted_peptides:
            if pep not in peptides_set:
                continue
            accepted_peptides.append(pep)
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
        logger.info(f"Accepted {len(accepted_peptides)} non-redundant peptides.")
        return accepted_peptides, redundant_mapping

    def _build_full_contig_map(self, peptides: List[str]) -> Dict[int, List[str]]:
        """Group layout and contained peptides by layout index."""
        full_contig_map = defaultdict(list)
        for pep in peptides:
            contig_idx = self.peptide_to_contig.get(pep, None)
            if contig_idx is not None:
                full_contig_map[contig_idx].append(pep)
        return full_contig_map

    def _calculate_overlap_contig_features(self, peptides: List[str]) -> pd.DataFrame:
        """Build layouts and calculate per-peptide layout features."""
        calculated_columns = [
            "clean_peptide",
            "contig_member_count",
            "contig_length",
            "log_contig_member_count",
            "contig_member_rank",
            "log_contig_member_rank",
            "contig_seq_length_diff",
            "contig_extension_ratio",
        ]
        if not peptides:
            self._overlap_graph = nx.DiGraph()
            self._simplified_graph = nx.DiGraph()
            self.assembled_contigs = []
            self.peptide_to_contig = {}
            return pd.DataFrame(columns=calculated_columns)

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

        full_contig_map = self._build_full_contig_map(peptides)

        for contig_index, full_peptides in full_contig_map.items():
            self.assembled_contigs[contig_index]["full_contig_peptides"] = full_peptides

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
        # TODO: Consider a true contig-level dense rank in a future version
        # Keep the peptide-row competition rank for v0.2.0 output compatibility
        features_df["contig_member_rank"] = features_df["contig_member_count"].rank(
            method="min",
            ascending=False,  # method may turn to "dense" in the future
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
        """Calculate features once and align them to the original peptide rows."""
        if self.overlap_data is None:
            self.overlap_data = pd.DataFrame(self.original_peptides, columns=["Peptide"])
            self.overlap_data["clean_peptide"] = self.overlap_data["Peptide"].apply(
                self._preprocess_peptides
            )
            self.filtered_peptides = self._filter_peptides(
                self.overlap_data["clean_peptide"].unique().tolist()
            )

            features_df = self._calculate_overlap_contig_features(self.filtered_peptides)

            logger.info("Mapping features back to original peptides.")
            self.overlap_data = self.overlap_data.merge(
                features_df, on="clean_peptide", how="left"
            )

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
        """Return public overlap-layout feature columns."""
        features_df = self._integrate_overlap_features()
        features_df = features_df[["Peptide"] + self.feature_columns]
        logger.info(f"Generated overlap features for {len(features_df)} peptides.")
        return features_df

    def get_full_data(self) -> pd.DataFrame:
        """Return legacy layout metadata, including contained peptide members."""
        self._integrate_overlap_features()
        if self.full_data is not None:
            logger.info("Full data has already been computed. Returning cached data.")
            return self.full_data
        data_list = []

        for peptide in tqdm(self.filtered_peptides):
            contig_idx = self.peptide_to_contig.get(peptide, None)
            if contig_idx is not None:
                contig_info = self.assembled_contigs[contig_idx]
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
            peptides=sorted(set(psms.df["sequence"])),
            min_overlap_length=params.get("minOverlapLength", 8),
            min_length=params.get("minLength", 8),
            max_length=params.get("maxLength", 25),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
        )
        instance._overlapping_score = params.get("overlappingScore", None)
        return instance

    def apply(self, psms):
        features = self.generate_features()
        full_data = self.get_full_data()
        features = features.rename(columns={"Peptide": "sequence"})
        psms.add_features(
            features,
            on="sequence",
            columns=self.feature_columns,
        )

        if self._overlapping_score:
            assign_brother_aggregated_feature(
                psms,
                full_data=full_data,
                feature_columns=self._overlapping_score,
            )

    def feature_groups(self, name):
        groups = super().feature_groups(name)
        if self._overlapping_score:
            groups["ContigFeatures"] = tuple(_contig_feature_columns(self._overlapping_score))
        return groups


def assign_brother_aggregated_feature(
    psms: PsmContainer,
    full_data: pd.DataFrame,
    feature_columns: Union[str, List[str]],
) -> None:
    """Add mean and count-scaled PSM features grouped by legacy layout sequence."""
    if isinstance(feature_columns, str):
        feature_columns = [feature_columns]
    psms_df = psms.df.merge(
        full_data[["Peptide", "ContigSequence"]].rename(columns={"Peptide": "sequence"}),
        on="sequence",
        how="left",
        validate="many_to_one",
    )

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

    agg_feature_columns = _contig_feature_columns(feature_columns)

    new_features_df = psms_with_agg[["psm_id", *agg_feature_columns]]
    psms.add_features(new_features_df, on="psm_id", columns=agg_feature_columns)


def _contig_feature_columns(feature_columns: Union[str, List[str]]) -> List[str]:
    if isinstance(feature_columns, str):
        feature_columns = [feature_columns]
    return [f"{feature}_contig_avg" for feature in feature_columns] + [
        f"{feature}_contig_sum" for feature in feature_columns
    ]


feature_generator_factory.register_generator(
    "OverlappingPeptide", OverlappingPeptideFeatureGenerator
)
