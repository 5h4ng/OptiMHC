# TODO: Except 'best' mode, the other modes seems to be not working properly. We need to investigate this issue.

import logging
from multiprocessing import Pool, cpu_count
from typing import List, Optional

import pandas as pd
from mhctools import NetMHCpan41
from tqdm import tqdm

from optimhc import utils
from optimhc.feature.factory import feature_generator_factory

from .base_feature_generator import BaseFeatureGenerator

logger = logging.getLogger(__name__)


# Each worker process gets its own copy of this global.
# Pool(initializer=_init_worker) sets it once per process,
# so netMHCpan is called O(processes) not O(chunks).
_worker_predictor: Optional[NetMHCpan41] = None


def _init_worker(alleles: List[str], program_name: str) -> None:
    global _worker_predictor
    _worker_predictor = NetMHCpan41(alleles=alleles, program_name=program_name)


def _predict_peptide_chunk(peptides_chunk: List[str]) -> pd.DataFrame:
    return _worker_predictor.predict_peptides(peptides_chunk).to_dataframe()


class NetMHCpanFeatureGenerator(BaseFeatureGenerator):
    """
    Generate NetMHCpan features for peptides based on specified MHC class I alleles.

    This generator calculates NetMHCpan binding predictions for each peptide against
    the provided MHC class I alleles.

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC allele names (e.g., ['HLA-A*02:01', 'HLA-B*07:02']).
    mode : str, optional
        Mode of feature generation. Options:
        - 'best': Return only the best allele information for each peptide.
        - 'all': Return predictions for all alleles with allele-specific suffixes plus best allele info.
        Default is 'best'.
    remove_pre_nxt_aa : bool, optional
        Whether to remove the previous and next amino acids from peptides. Default is False.
    remove_modification : bool, optional
        Whether to remove modifications from peptides. Default is True.
    n_processes : int, optional
        Number of processes to use for multiprocessing. Default is 1.
    show_progress : bool, optional
        Whether to display a progress bar. Default is False.
    executable_path : str, optional
        Path to the netMHCpan executable. Defaults to "netMHCpan" (PATH lookup).

    Notes
    -----
    The generated features include:
    - netmhcpan_score: Raw binding score
    - netmhcpan_affinity: Binding affinity in nM
    - netmhcpan_percentile_rank: Percentile rank of the binding score
    """

    MIN_PEPTIDE_LENGTH = 8
    MAX_PEPTIDE_LENGTH = 30
    CHUNKSIZE = 250

    def __init__(
        self,
        peptides: List[str],
        alleles: List[str],
        mode: str = "best",
        remove_pre_nxt_aa: bool = False,
        remove_modification: bool = True,
        n_processes: int = 1,
        show_progress: bool = False,
        executable_path: str = "netMHCpan",
    ):
        if mode not in ["best", "all"]:
            raise ValueError("Mode must be one of 'best' or 'all'.")

        self.peptides = peptides
        self.alleles = alleles
        self.mode = mode
        if len(alleles) == 1:
            self.mode = "best"
            logger.info("Only one allele provided. Switching to 'best' mode.")
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.remove_modification = remove_modification
        self.n_processes = min(n_processes, cpu_count())
        self.show_progress = show_progress
        self.executable_path = executable_path
        self.predictor = (
            NetMHCpan41(alleles=self.alleles, program_name=executable_path)
            if self.n_processes == 1
            else None
        )
        self.predictions = None
        self._raw_predictions = None
        logger.info(
            f"Initialized NetMHCpanFeatureGenerator with {len(peptides)} peptides, alleles: {alleles}, mode: {mode}, n_processes: {self.n_processes}, show_progress: {self.show_progress}"
        )

    @property
    def feature_columns(self) -> List[str]:
        """
        Return the list of generated feature column names, determined by the mode.
        Only includes numerical features, excluding any string features like allele names.

        Returns
        -------
        List[str]
            List of feature column names:
            - For 'all' mode: netmhcpan_score_{allele}, netmhcpan_affinity_{allele},
              netmhcpan_percentile_rank_{allele} for each allele
            - For both modes: netmhcpan_best_score, netmhcpan_best_affinity,
              netmhcpan_best_percentile_rank
        """
        columns = []
        if self.mode == "all":
            allele_specific = []
            for allele in self.alleles:
                allele_specific.extend(
                    [
                        f"netmhcpan_score_{allele}",
                        f"netmhcpan_affinity_{allele}",
                        f"netmhcpan_percentile_rank_{allele}",
                    ]
                )
            columns.extend(allele_specific)

        # Both 'best' and 'all' modes include best allele numerical information
        columns.extend(
            [
                "netmhcpan_best_score",
                "netmhcpan_best_affinity",
                "netmhcpan_best_percentile_rank",
            ]
        )
        return columns

    @property
    def id_column(self) -> List[str]:
        return ["Peptide"]

    def _preprocess_peptides(self, peptide: str) -> str:
        if self.remove_pre_nxt_aa:
            peptide = utils.strip_flanking_and_charge(peptide)
        if self.remove_modification:
            peptide = utils.remove_modifications(peptide)
        peptide = peptide.replace("U", "C")
        return peptide

    def _predict_multiprocessing(self, peptides_to_predict: List[str]) -> pd.DataFrame:
        """
        Run NetMHCpan predictions using multiprocessing.

        One predictor is initialized per worker process (via Pool initializer),
        so netMHCpan -listMHC is called once per process rather than once per chunk.
        """
        logger.info("Running NetMHCpan predictions with multiprocessing.")
        chunk_size = min(
            NetMHCpanFeatureGenerator.CHUNKSIZE,
            max(1, len(peptides_to_predict) // self.n_processes),
        )
        chunks = [
            peptides_to_predict[i : i + chunk_size]
            for i in range(0, len(peptides_to_predict), chunk_size)
        ]

        with Pool(
            processes=self.n_processes,
            initializer=_init_worker,
            initargs=(self.alleles, self.executable_path),
        ) as pool:
            it = pool.imap(_predict_peptide_chunk, chunks)
            if self.show_progress:
                it = tqdm(it, total=len(chunks), desc="Predicting NetMHCpan")
            results = list(it)

        logger.info(
            f"Completed multiprocessing predictions for {len(peptides_to_predict)} peptides."
        )
        return pd.concat(results, ignore_index=True)

    def _predict(self) -> pd.DataFrame:
        if self.predictions is not None:
            logger.info("NetMHCpan predictions already exist. Skipping prediction.")
            return self.predictions

        logger.info("Starting NetMHCpan predictions.")
        predictions = pd.DataFrame(self.peptides, columns=["Peptide"])
        predictions["clean_peptide"] = predictions["Peptide"].apply(self._preprocess_peptides)

        peptides_to_predict = (
            predictions[
                predictions["clean_peptide"].apply(
                    lambda x: (
                        NetMHCpanFeatureGenerator.MIN_PEPTIDE_LENGTH
                        <= len(x)
                        <= NetMHCpanFeatureGenerator.MAX_PEPTIDE_LENGTH
                    )
                )
            ]["clean_peptide"]
            .unique()
            .tolist()
        )

        logger.info(
            f"Found {len(peptides_to_predict)} unique peptides meeting the length requirements."
        )

        if self.n_processes > 1:
            netmhcpan_results = self._predict_multiprocessing(peptides_to_predict)
        else:
            netmhcpan_results = self.predictor.predict_peptides(peptides_to_predict).to_dataframe()

        self._raw_predictions = netmhcpan_results.copy()
        logger.info(f"Predicted NetMHCpan results for {len(netmhcpan_results)} peptides.")

        predictions = predictions.merge(
            netmhcpan_results, left_on="clean_peptide", right_on="peptide", how="left"
        )
        predictions.drop(columns=["clean_peptide"], inplace=True)

        self.predictions = predictions
        logger.info(f"Completed NetMHCpan predictions for {len(peptides_to_predict)} peptides.")
        return self.predictions

    @property
    def raw_predictions(self) -> pd.DataFrame:
        if self._raw_predictions is None:
            self._predict()
        return self._raw_predictions

    def generate_features(self) -> pd.DataFrame:
        """
        Generate the final feature table with NetMHCpan features for each peptide.

        Returns
        -------
        pd.DataFrame
            DataFrame containing peptides and their predicted features.

        Notes
        -----
        The features generated depend on the mode:
        - 'best': Only the best allele information for each peptide
        - 'all': All allele predictions plus best allele information

        Missing values are filled with column medians.
        """
        predictions_df = self._predict()

        features_df = pd.DataFrame({"Peptide": self.peptides})

        if self.mode == "all":
            features_df = self._generate_all_allele_features(predictions_df, features_df)
        features_df = self._generate_best_allele_features(predictions_df, features_df)

        features_df = self._fill_missing_values(features_df)

        selected_columns = ["Peptide"] + self.feature_columns
        logger.info(f"Final selected feature columns: {selected_columns}")
        features_df = features_df[selected_columns]

        if features_df.isna().sum().sum() > 0:
            logger.warning(
                "NaN values still exist in the generated features after filling with median/mode values."
            )

        return features_df

    def _generate_all_allele_features(
        self, predictions_df: pd.DataFrame, features_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Generate features for all alleles.

        Parameters
        ----------
        predictions_df : pd.DataFrame
            The predictions DataFrame.
        features_df : pd.DataFrame
            The features DataFrame to update.

        Returns
        -------
        pd.DataFrame
            Updated features DataFrame with all allele features:
            - Peptide: Original peptide sequence
            - netmhcpan_score_{allele}: Raw binding score for each allele
            - netmhcpan_affinity_{allele}: Binding affinity for each allele
            - netmhcpan_percentile_rank_{allele}: Percentile rank for each allele
        """
        logger.info("Generating features for all alleles.")

        for allele in self.alleles:
            logger.info(f"Adding scores for allele {allele}.")
            allele_df = predictions_df[predictions_df["allele"] == allele].copy()

            if allele_df.empty:
                logger.warning(
                    f"No prediction results found for allele {allele}. Filling with NaN."
                )
                allele_features = pd.DataFrame(
                    {
                        "Peptide": self.peptides,
                        f"netmhcpan_score_{allele}": [pd.NA] * len(self.peptides),
                        f"netmhcpan_affinity_{allele}": [pd.NA] * len(self.peptides),
                        f"netmhcpan_percentile_rank_{allele}": [pd.NA] * len(self.peptides),
                    }
                )
            else:
                allele_df = allele_df.rename(
                    columns={
                        "score": f"netmhcpan_score_{allele}",
                        "affinity": f"netmhcpan_affinity_{allele}",
                        "percentile_rank": f"netmhcpan_percentile_rank_{allele}",
                    }
                )
                allele_features = allele_df[
                    [
                        "Peptide",
                        f"netmhcpan_score_{allele}",
                        f"netmhcpan_affinity_{allele}",
                        f"netmhcpan_percentile_rank_{allele}",
                    ]
                ]

            features_df = features_df.merge(allele_features, on="Peptide", how="left")

        logger.info("Added scores for all alleles.")
        return features_df

    def _generate_best_allele_features(
        self, predictions_df: pd.DataFrame, features_df: pd.DataFrame
    ) -> pd.DataFrame:
        logger.info("Generating features for best allele.")

        valid_predictions = predictions_df.dropna(subset=["percentile_rank"])

        if valid_predictions.empty:
            logger.warning("No valid predictions available to select the best allele.")
            best_allele_features = pd.DataFrame(
                {
                    "Peptide": self.peptides,
                    "netmhcpan_best_allele": ["Unknown"] * len(self.peptides),
                    "netmhcpan_best_score": [pd.NA] * len(self.peptides),
                    "netmhcpan_best_affinity": [pd.NA] * len(self.peptides),
                    "netmhcpan_best_percentile_rank": [pd.NA] * len(self.peptides),
                }
            )
        else:
            idx = valid_predictions.groupby("Peptide")["percentile_rank"].idxmin()

            best_allele_features = valid_predictions.loc[idx].rename(
                columns={
                    "allele": "netmhcpan_best_allele",
                    "score": "netmhcpan_best_score",
                    "affinity": "netmhcpan_best_affinity",
                    "percentile_rank": "netmhcpan_best_percentile_rank",
                }
            )
            best_allele_features = best_allele_features[
                [
                    "Peptide",
                    "netmhcpan_best_allele",
                    "netmhcpan_best_score",
                    "netmhcpan_best_affinity",
                    "netmhcpan_best_percentile_rank",
                ]
            ]

            missing_peptides = set(self.peptides) - set(best_allele_features["Peptide"])
            if missing_peptides:
                logger.warning(
                    f"Found {len(missing_peptides)} peptides with no best allele prediction."
                )
                missing_features = pd.DataFrame(
                    {
                        "Peptide": list(missing_peptides),
                        "netmhcpan_best_allele": ["Unknown"] * len(missing_peptides),
                        "netmhcpan_best_score": [pd.NA] * len(missing_peptides),
                        "netmhcpan_best_affinity": [pd.NA] * len(missing_peptides),
                        "netmhcpan_best_percentile_rank": [pd.NA] * len(missing_peptides),
                    }
                )
                best_allele_features = pd.concat(
                    [best_allele_features, missing_features], ignore_index=True
                )

        features_df = features_df.merge(best_allele_features, on="Peptide", how="left")
        logger.info("Added best allele information.")
        return features_df

    def _fill_missing_values(self, features_df: pd.DataFrame) -> pd.DataFrame:
        logger.info("Filling missing values in the features DataFrame.")

        if "netmhcpan_best_allele" in features_df.columns:
            features_df["netmhcpan_best_allele"] = features_df["netmhcpan_best_allele"].fillna(
                "Unknown"
            )

        if self.mode == "all":
            for allele in self.alleles:
                for metric in ["score", "affinity", "percentile_rank"]:
                    col = f"netmhcpan_{metric}_{allele}"
                    if col in features_df.columns:
                        features_df[col] = features_df[col].fillna(features_df[col].median())

        for metric in ["best_score", "best_affinity", "best_percentile_rank"]:
            col = f"netmhcpan_{metric}"
            if col in features_df.columns and features_df[col].isna().any():
                median_value = features_df[col].median()
                features_df[col] = features_df[col].fillna(
                    0 if pd.isna(median_value) else median_value
                )

        logger.info("Filled missing values in the features DataFrame.")
        return features_df

    @classmethod
    def from_config(cls, psms, config, params):
        return cls(
            peptides=list(set(psms.peptides)),
            alleles=config.get("allele", []),
            mode=params.get("mode", "best"),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
            n_processes=config.get("numProcesses", 1),
            show_progress=config.get("showProgress", False),
            executable_path=params.get("executablePath", "netMHCpan"),
        )


feature_generator_factory.register_generator("NetMHCpan", NetMHCpanFeatureGenerator)
