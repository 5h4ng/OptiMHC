# TODO: Except 'best' mode, the other modes seems to be not working properly. We need to investigate this issue.

import logging
from functools import partial
from multiprocessing import Pool, cpu_count
from typing import List

import pandas as pd
from mhctools import NetMHCpan41
from tqdm import tqdm

from optimhc import utils

from .base_feature_generator import BaseFeatureGenerator

logger = logging.getLogger(__name__)


# Helper function for multiprocessing
def _predict_peptide_chunk(peptides_chunk: List[str], alleles: List[str]) -> pd.DataFrame:
    """
    Predict NetMHCpan scores for a chunk of peptides.

    Parameters
    ----------
    peptides_chunk : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC allele names.

    Returns
    -------
    pd.DataFrame
        DataFrame containing predictions:
        - peptide: Peptide sequence
        - allele: MHC allele
        - score: Raw binding score
        - affinity: Binding affinity in nM
        - percentile_rank: Percentile rank
    """
    predictor = NetMHCpan41(alleles=alleles)
    results = predictor.predict_peptides(peptides_chunk)
    return results.to_dataframe()


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
        Whether to include the previous and next amino acids in peptides.
        If True, remove them. Default is True.
    remove_modification : bool, optional
        Whether to include modifications in peptides.
        If True, remove them. Default is True.
    n_processes : int, optional
        Number of processes to use for multiprocessing.
        Default is 1 (no multiprocessing).
    show_progress : bool, optional
        Whether to display a progress bar. Default is False.

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
        self.predictor = NetMHCpan41(alleles=self.alleles)
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
        Run NetMHCpan predictions using multiprocessing.

        Parameters
        ----------
        peptides_to_predict : List[str]
            List of peptides to predict.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the prediction results:
            - peptide: Peptide sequence
            - allele: MHC allele
            - score: Raw binding score
            - affinity: Binding affinity in nM
            - percentile_rank: Percentile rank

        Notes
        -----
        This method:
        1. Splits peptides into chunks
        2. Processes chunks in parallel
        3. Combines results into a single DataFrame
        """
        logger.info("Running NetMHCpan predictions with multiprocessing.")
        chunksize = min(
            NetMHCpanFeatureGenerator.CHUNKSIZE,
            max(1, len(peptides_to_predict) // self.n_processes),
        )
        func = partial(_predict_peptide_chunk, alleles=self.alleles)

        with Pool(processes=self.n_processes) as pool:
            if self.show_progress:
                results = list(
                    tqdm(
                        pool.imap(
                            func,
                            [
                                peptides_to_predict[i : i + chunksize]
                                for i in range(0, len(peptides_to_predict), chunksize)
                            ],
                        ),
                        total=(len(peptides_to_predict) + chunksize - 1) // chunksize,
                        desc="Predicting NetMHCpan",
                    )
                )
            else:
                results = pool.map(
                    func,
                    [
                        peptides_to_predict[i : i + chunksize]
                        for i in range(0, len(peptides_to_predict), chunksize)
                    ],
                )

        netmhcpan_results = pd.concat(results, ignore_index=True)
        # Save the raw prediction results
        self._raw_predictions = netmhcpan_results.copy()
        logger.info(
            f"Completed multiprocessing predictions for {len(peptides_to_predict)} peptides."
        )
        return netmhcpan_results

    def _predict(self) -> pd.DataFrame:
        """
        Run NetMHCpan predictions and cache the result.

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
        2. Filters peptides by length (8-30 amino acids)
        3. Runs predictions (with or without multiprocessing)
        4. Merges results with original peptides
        """
        if self.predictions is not None:
            logger.info("NetMHCpan predictions already exist. Skipping prediction.")
            return self.predictions

        logger.info("Starting NetMHCpan predictions.")
        self.predictions = pd.DataFrame(self.peptides, columns=["Peptide"])
        self.predictions["clean_peptide"] = self.predictions["Peptide"].apply(
            self._preprocess_peptides
        )

        # Filter peptides that meet the length requirements
        peptides_to_predict = (
            self.predictions[
                self.predictions["clean_peptide"].apply(
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
            # If not using multiprocessing, save raw prediction results here
            self._raw_predictions = netmhcpan_results.copy()

        logger.info(f"Predicted NetMHCpan results for {len(netmhcpan_results)} peptides.")

        self.predictions = self.predictions.merge(
            netmhcpan_results, left_on="clean_peptide", right_on="peptide", how="left"
        )
        self.predictions.drop(columns=["clean_peptide"], inplace=True)

        logger.info(f"Completed NetMHCpan predictions for {len(peptides_to_predict)} peptides.")
        return self.predictions

    @property
    def raw_predictions(self) -> pd.DataFrame:
        """
        Return the raw prediction results from NetMHCpan.

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
        Get the raw prediction results DataFrame from NetMHCpan.

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

    def generate_features(self) -> pd.DataFrame:
        """
        Generate the final feature table with NetMHCpan features for each peptide.

        Returns
        -------
        pd.DataFrame
            DataFrame containing peptides and their predicted features:
            - Peptide: Original peptide sequence
            - For 'all' mode: netmhcpan_score_{allele}, netmhcpan_affinity_{allele},
              netmhcpan_percentile_rank_{allele} for each allele
            - For both modes: netmhcpan_best_score, netmhcpan_best_affinity,
              netmhcpan_best_percentile_rank

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
            - netmhcpan_best_allele: Best binding allele
            - netmhcpan_best_score: Best binding score
            - netmhcpan_best_affinity: Best binding affinity
            - netmhcpan_best_percentile_rank: Best percentile rank

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
                    "netmhcpan_best_allele": ["Unknown"] * len(self.peptides),
                    "netmhcpan_best_score": [pd.NA] * len(self.peptides),
                    "netmhcpan_best_affinity": [pd.NA] * len(self.peptides),
                    "netmhcpan_best_percentile_rank": [pd.NA] * len(self.peptides),
                }
            )
        else:
            # Find the index of the minimum percentile rank for each peptide
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

            # Handle missing peptides
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
        if "netmhcpan_best_allele" in features_df.columns:
            features_df["netmhcpan_best_allele"].fillna("Unknown", inplace=True)

        # Fill numeric values with median for all allele features
        if self.mode == "all":
            for allele in self.alleles:
                for metric in ["score", "affinity", "percentile_rank"]:
                    col = f"netmhcpan_{metric}_{allele}"
                    if col in features_df.columns:
                        median_value = features_df[col].median()
                        features_df[col].fillna(median_value, inplace=True)

        # Fill numeric values for best allele features
        for metric in ["best_score", "best_affinity", "best_percentile_rank"]:
            col = f"netmhcpan_{metric}"
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

    # def get_best_allele(self) -> pd.DataFrame:
    #     """
    #     Return the best allele (with the lowest percentile rank) for each peptide across all alleles.

    #     Returns:
    #         pd.DataFrame: DataFrame containing the best allele information.
    #     """
    #     logger.info("Getting best allele information.")
    #     predictions_df = self._predict()

    #     features_df = pd.DataFrame({'Peptide': self.peptides})
    #     best_features_df = self._generate_best_allele_features(predictions_df, features_df)
    #     best_columns = ['Peptide', 'netmhcpan_best_allele', 'netmhcpan_best_score',
    #                      'netmhcpan_best_affinity', 'netmhcpan_best_percentile_rank']
    #     best_allele_df = best_features_df[best_columns]
    #     best_allele_df = self._fill_missing_values(best_allele_df)

    #     logger.info(f"Generated best allele information for {len(best_allele_df)} peptides.")
    #     return best_allele_df
