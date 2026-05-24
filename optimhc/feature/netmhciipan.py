# TODO: set 'BA' and 'EL' as optional parameters for the user to choose the prediction method.

import logging
import os
from multiprocessing import Pool, cpu_count
from typing import List, Optional

import pandas as pd
from mhctools import NetMHCIIpan43_BA
from tqdm import tqdm

from optimhc import utils
from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)


# Each worker process gets its own copy of this global.
# Pool(initializer=_init_worker) sets it once per process,
# so netMHCIIpan is called O(processes) not O(chunks).
_worker_predictor: Optional[NetMHCIIpan43_BA] = None


def _init_worker(alleles: List[str], program_name: str) -> None:
    global _worker_predictor
    # mhctools uses program_name as a tempfile prefix — pass basename, add dir to PATH.
    exe_dir = os.path.dirname(program_name)
    if exe_dir:
        os.environ["PATH"] = exe_dir + os.pathsep + os.environ.get("PATH", "")
        program_name = os.path.basename(program_name)
    _worker_predictor = NetMHCIIpan43_BA(alleles=alleles, program_name=program_name)


def _predict_peptide_chunk(peptides_chunk: List[str]) -> pd.DataFrame:
    return _worker_predictor.predict_peptides(peptides_chunk).to_dataframe()


class NetMHCIIpanFeatureGenerator(BaseFeatureGenerator):
    """
    Generate NetMHCIIpan features for given peptides based on specified MHC Class II alleles.

    This feature generator uses the NetMHCIIpan43_BA interface to predict MHC Class II binding
    for each peptide and returns scores and features based on the specified parameters.

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC Class II alleles, e.g., ['DRB1_0101', 'DRB1_0102'].
    mode : str, optional
        Feature generation mode. Options:
        - 'best': Return only the best result for each peptide across all alleles.
        - 'all': Return prediction results for each peptide across all alleles (with allele-specific column suffixes).
        Default is 'best'.
    remove_pre_nxt_aa : bool, optional
        Whether to remove the amino acids flanking the peptide (e.g., removing X-AA/AA-X forms).
        Default is True.
    remove_modification : bool, optional
        Whether to remove modification information from peptides, e.g., (Phospho).
        Default is True.
    n_processes : int, optional
        Number of processes to use. Default is 1 (no multiprocessing).
    show_progress : bool, optional
        Whether to display a progress bar. Default is False.
    executable_path : str, optional
        Path to the netMHCIIpan executable. Defaults to "netMHCIIpan" (PATH lookup).

    Notes
    -----
    The generated features include:
    - netmhciipan_score: Raw binding score
    - netmhciipan_affinity: Binding affinity in nM
    - netmhciipan_percentile_rank: Percentile rank of the binding score
    """

    MIN_PEPTIDE_LENGTH = 9  # Minimum peptide length for MHC Class II is usually 9
    MAX_PEPTIDE_LENGTH = 50  # Can be adjusted based on use case
    CHUNKSIZE = 500

    def __init__(
        self,
        peptides: List[str],
        alleles: List[str],
        mode: str = "best",
        remove_pre_nxt_aa: bool = True,
        remove_modification: bool = True,
        n_processes: int = 1,
        show_progress: bool = False,
        executable_path: str = "netMHCIIpan",
        *args,
        **kwargs,
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
        # mhctools uses program_name as a tempfile prefix — pass basename, add dir to PATH.
        exe_dir = os.path.dirname(executable_path)
        if exe_dir:
            os.environ["PATH"] = exe_dir + os.pathsep + os.environ.get("PATH", "")
            executable_path = os.path.basename(executable_path)
        self.predictor = NetMHCIIpan43_BA(alleles=self.alleles, program_name=executable_path)
        self.predictions = None
        self._raw_predictions = None

        logger.info(
            f"Initialized NetMHCIIpanFeatureGenerator with {len(peptides)} peptides, "
            f"alleles={alleles}, mode={self.mode}, "
            f"n_processes={self.n_processes}, show_progress={self.show_progress}"
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
            - For 'all' mode: netmhciipan_score_{allele}, netmhciipan_affinity_{allele},
              netmhciipan_percentile_rank_{allele} for each allele
            - For both modes: netmhciipan_best_score, netmhciipan_best_affinity,
              netmhciipan_best_percentile_rank
        """
        columns = []
        if self.mode == "all":
            allele_specific = []
            for allele in self.alleles:
                allele_specific.extend(
                    [
                        f"netmhciipan_score_{allele}",
                        f"netmhciipan_affinity_{allele}",
                        f"netmhciipan_percentile_rank_{allele}",
                    ]
                )
            columns.extend(allele_specific)

        # Both 'best' and 'all' modes include best allele numerical information
        columns.extend(
            [
                "netmhciipan_best_score",
                "netmhciipan_best_affinity",
                "netmhciipan_best_percentile_rank",
            ]
        )
        return columns

    @property
    def id_column(self) -> List[str]:
        """
        Return the list of input columns required for the feature generator.

        Returns
        -------
        List[str]
            List of input column names.
        """
        return ["Peptide"]

    def _preprocess_peptides(self, peptide: str) -> str:
        """
        Preprocess the input peptide by removing flanking amino acids,
        modifications, and replacing non-standard amino acids.

        Parameters
        ----------
        peptide : str
            The original peptide sequence.

        Returns
        -------
        str
            The preprocessed peptide sequence.
        """
        if self.remove_pre_nxt_aa:
            peptide = utils.strip_flanking_and_charge(peptide)
        if self.remove_modification:
            peptide = utils.remove_modifications(peptide)
        # Replace any non-standard amino acids if necessary
        peptide = peptide.replace("U", "C")
        return peptide

    def _predict_multiprocessing(self, peptides_to_predict: List[str]) -> pd.DataFrame:
        """
        Run NetMHCIIpan predictions using multiprocessing.

        One predictor is initialized per worker process (via Pool initializer),
        so netMHCIIpan -list is called once per process rather than once per chunk.

        Parameters
        ----------
        peptides_to_predict : List[str]
            List of peptides to predict.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns: peptide, allele, score, affinity, percentile_rank.
        """
        logger.info("Running NetMHCIIpan predictions with multiprocessing.")
        chunk_size = min(
            NetMHCIIpanFeatureGenerator.CHUNKSIZE,
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
                it = tqdm(it, total=len(chunks), desc="Predicting NetMHCIIpan")
            results = list(it)

        logger.info(
            f"Completed multiprocessing predictions for {len(peptides_to_predict)} peptides."
        )
        return pd.concat(results, ignore_index=True)

    def _predict(self) -> pd.DataFrame:
        """
        Run NetMHCIIpan predictions and cache the result.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the prediction results:
            - Peptide: Original peptide sequence
            - peptide: Cleaned peptide sequence
            - allele: MHC allele
            - score: Raw binding score
            - affinity: Binding affinity in nM
            - percentile_rank: Percentile rank

        Notes
        -----
        This method:
        1. Preprocesses peptides
        2. Filters peptides by length (9-50 amino acids)
        3. Runs predictions (with or without multiprocessing)
        4. Merges results with original peptides
        """
        if self.predictions is not None:
            logger.info("NetMHCIIpan predictions already exist. Skipping prediction.")
            return self.predictions

        logger.info("Starting NetMHCIIpan predictions.")
        self.predictions = pd.DataFrame(self.peptides, columns=["Peptide"])
        self.predictions["clean_peptide"] = self.predictions["Peptide"].apply(
            self._preprocess_peptides
        )

        # Filter peptides that meet the length requirements
        peptides_to_predict = (
            self.predictions[
                self.predictions["clean_peptide"].apply(
                    lambda x: (
                        NetMHCIIpanFeatureGenerator.MIN_PEPTIDE_LENGTH
                        <= len(x)
                        <= NetMHCIIpanFeatureGenerator.MAX_PEPTIDE_LENGTH
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
            netmhciipan_results = self._predict_multiprocessing(peptides_to_predict)
        else:
            netmhciipan_results = self.predictor.predict_peptides(
                peptides_to_predict
            ).to_dataframe()

        self._raw_predictions = netmhciipan_results.copy()
        logger.info(f"Predicted NetMHCIIpan results for {len(netmhciipan_results)} peptides.")

        self.predictions = self.predictions.merge(
            netmhciipan_results, left_on="clean_peptide", right_on="peptide", how="left"
        )
        self.predictions.drop(columns=["clean_peptide"], inplace=True)

        logger.info(f"Completed NetMHCIIpan predictions for {len(peptides_to_predict)} peptides.")
        return self.predictions

    def generate_features(self) -> pd.DataFrame:
        """
        Generate the final feature table with NetMHCIIpan features for each peptide.

        Returns
        -------
        pd.DataFrame
            DataFrame containing peptides and their predicted features:
            - Peptide: Original peptide sequence
            - For 'all' mode: netmhciipan_score_{allele}, netmhciipan_affinity_{allele},
              netmhciipan_percentile_rank_{allele} for each allele
            - For both modes: netmhciipan_best_score, netmhciipan_best_affinity,
              netmhciipan_best_percentile_rank

        Notes
        -----
        The features generated depend on the mode:
        - 'best': Only the best allele information for each peptide
        - 'all': All allele predictions plus best allele information

        Missing values are handled consistently by filling with median values
        for numeric columns.
        """
        predictions_df = self._predict()

        features_df = pd.DataFrame({"Peptide": self.peptides})

        # Generate allele-specific features if mode is 'all', otherwise generate best allele features
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
            - netmhciipan_score_{allele}: Raw binding score for each allele
            - netmhciipan_affinity_{allele}: Binding affinity for each allele
            - netmhciipan_percentile_rank_{allele}: Percentile rank for each allele
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
                        f"netmhciipan_score_{allele}": [pd.NA] * len(self.peptides),
                        f"netmhciipan_affinity_{allele}": [pd.NA] * len(self.peptides),
                        f"netmhciipan_percentile_rank_{allele}": [pd.NA] * len(self.peptides),
                    }
                )
            else:
                allele_df = allele_df.rename(
                    columns={
                        "score": f"netmhciipan_score_{allele}",
                        "affinity": f"netmhciipan_affinity_{allele}",
                        "percentile_rank": f"netmhciipan_percentile_rank_{allele}",
                    }
                )

                allele_features = allele_df[
                    [
                        "Peptide",
                        f"netmhciipan_score_{allele}",
                        f"netmhciipan_affinity_{allele}",
                        f"netmhciipan_percentile_rank_{allele}",
                    ]
                ]

            features_df = features_df.merge(allele_features, on="Peptide", how="left")

        logger.info("Added scores for all alleles.")
        return features_df

    def _generate_best_allele_features(
        self, predictions_df: pd.DataFrame, features_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Generate features for the best allele.

        Parameters
        ----------
        predictions_df : pd.DataFrame
            The predictions DataFrame.
        features_df : pd.DataFrame
            The features DataFrame to update.

        Returns
        -------
        pd.DataFrame
            Updated features DataFrame with best allele features:
            - Peptide: Original peptide sequence
            - netmhciipan_best_allele: Best binding allele
            - netmhciipan_best_score: Best binding score
            - netmhciipan_best_affinity: Best binding affinity
            - netmhciipan_best_percentile_rank: Best percentile rank

        Notes
        -----
        The best allele is determined by the lowest percentile rank.
        """
        logger.info("Generating features for best allele.")

        valid_predictions = predictions_df.dropna(subset=["percentile_rank"])

        if valid_predictions.empty:
            logger.warning("No valid predictions available to select the best allele.")
            best_allele_features = pd.DataFrame(
                {
                    "Peptide": self.peptides,
                    "netmhciipan_best_allele": ["Unknown"] * len(self.peptides),
                    "netmhciipan_best_score": [pd.NA] * len(self.peptides),
                    "netmhciipan_best_affinity": [pd.NA] * len(self.peptides),
                    "netmhciipan_best_percentile_rank": [pd.NA] * len(self.peptides),
                }
            )
        else:
            # Find the index of the minimum percentile rank for each peptide
            idx = valid_predictions.groupby("Peptide")["percentile_rank"].idxmin()

            best_allele_features = valid_predictions.loc[idx].rename(
                columns={
                    "allele": "netmhciipan_best_allele",
                    "score": "netmhciipan_best_score",
                    "affinity": "netmhciipan_best_affinity",
                    "percentile_rank": "netmhciipan_best_percentile_rank",
                }
            )

            best_allele_features = best_allele_features[
                [
                    "Peptide",
                    "netmhciipan_best_allele",
                    "netmhciipan_best_score",
                    "netmhciipan_best_affinity",
                    "netmhciipan_best_percentile_rank",
                ]
            ]

            # Handle missing peptides
            missing_peptides = set(self.peptides) - set(best_allele_features["Peptide"])

            if missing_peptides:
                logger.warning(
                    f"Found {len(missing_peptides)} peptides with no best allele prediction."
                )
                missing_features = pd.DataFrame(
                    {
                        "Peptide": list(missing_peptides),
                        "netmhciipan_best_allele": ["Unknown"] * len(missing_peptides),
                        "netmhciipan_best_score": [pd.NA] * len(missing_peptides),
                        "netmhciipan_best_affinity": [pd.NA] * len(missing_peptides),
                        "netmhciipan_best_percentile_rank": [pd.NA] * len(missing_peptides),
                    }
                )
                best_allele_features = pd.concat(
                    [best_allele_features, missing_features], ignore_index=True
                )

        features_df = features_df.merge(best_allele_features, on="Peptide", how="left")
        logger.info("Added best allele information.")

        return features_df

    def _fill_missing_values(self, features_df: pd.DataFrame) -> pd.DataFrame:
        """
        Fill missing values in the features DataFrame.

        Parameters
        ----------
        features_df : pd.DataFrame
            The features DataFrame to update.

        Returns
        -------
        pd.DataFrame
            Updated features DataFrame with filled missing values.

        Notes
        -----
        This method:
        1. Fills best allele string values with 'Unknown'
        2. Fills numeric values with median for all allele features
        3. Fills numeric values for best allele features with median
        """
        logger.info("Filling missing values in the features DataFrame.")

        # Fill best allele string values
        if "netmhciipan_best_allele" in features_df.columns:
            features_df["netmhciipan_best_allele"].fillna("Unknown", inplace=True)

        # Fill numeric values with median for all allele features
        if self.mode == "all":
            for allele in self.alleles:
                for metric in ["score", "affinity", "percentile_rank"]:
                    col = f"netmhciipan_{metric}_{allele}"
                    if col in features_df.columns:
                        median_value = features_df[col].median()
                        features_df[col].fillna(median_value, inplace=True)

        # Fill numeric values for best allele features
        for metric in ["best_score", "best_affinity", "best_percentile_rank"]:
            col = f"netmhciipan_{metric}"
            if col in features_df.columns and features_df[col].isna().any():
                median_value = features_df[col].median()
                # If all values are NA, median will be NA, so use 0 instead
                median_value = 0 if pd.isna(median_value) else median_value
                features_df[col].fillna(median_value, inplace=True)

        logger.info("Filled missing values in the features DataFrame.")
        return features_df

    def predictions_to_dataframe(self) -> pd.DataFrame:
        """
        Convert the predictions to a DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the predictions.

        Raises
        ------
        ValueError
            If no predictions are available.
        """
        if self.predictions is None:
            raise ValueError("No predictions available. Please run 'generate_features' first.")
        return self.predictions

    @property
    def raw_predictions(self) -> pd.DataFrame:
        """
        Return the raw prediction results from NetMHCIIpan.

        Returns
        -------
        pd.DataFrame
            Raw prediction results DataFrame containing:
            - peptide: Cleaned peptide sequence
            - allele: MHC allele
            - score: Raw binding score
            - affinity: Binding affinity in nM
            - percentile_rank: Percentile rank
        """
        if self._raw_predictions is None:
            self._predict()
        return self._raw_predictions

    def get_raw_predictions(self) -> pd.DataFrame:
        """
        Get the raw prediction results DataFrame from NetMHCIIpan.

        Returns
        -------
        pd.DataFrame
            Raw prediction results DataFrame containing:
            - peptide: Cleaned peptide sequence
            - allele: MHC allele
            - score: Raw binding score
            - affinity: Binding affinity in nM
            - percentile_rank: Percentile rank
        """
        return self.raw_predictions

    def save_raw_predictions(self, file_path: str, **kwargs) -> None:
        """
        Save the raw prediction results to a file.

        Parameters
        ----------
        file_path : str
            Path to save the file.
        **kwargs : dict
            Additional parameters passed to pandas.DataFrame.to_csv.
            If 'index' is not specified, it defaults to False.

        Notes
        -----
        This method saves the raw predictions DataFrame to a CSV file.
        The DataFrame includes:
        - peptide: Cleaned peptide sequence
        - allele: MHC allele
        - score: Raw binding score
        - affinity: Binding affinity in nM
        - percentile_rank: Percentile rank
        """
        if "index" not in kwargs:
            kwargs["index"] = False
        if self.raw_predictions is not None:
            self.raw_predictions.to_csv(file_path, **kwargs)
            logger.info(f"Raw prediction results saved to: {file_path}")
        else:
            logger.warning("No raw prediction results available to save.")

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
            executable_path=params.get("executablePath", "netMHCIIpan"),
        )


feature_generator_factory.register_generator("NetMHCIIpan", NetMHCIIpanFeatureGenerator)
