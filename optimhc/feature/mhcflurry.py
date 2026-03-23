import logging
from typing import List

import pandas as pd
from mhcflurry import Class1PresentationPredictor

from optimhc import utils
from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)


class MHCflurryFeatureGenerator(BaseFeatureGenerator):
    """
    Generate MHCflurry features for peptides based on specified MHC class I alleles.

    This generator calculates MHCflurry presentation scores for each peptide against
    the provided MHC class I alleles.

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC allele names (e.g., ['HLA-A01:01', 'HLA-B07:02']).
    remove_pre_nxt_aa : bool, optional
        Whether to include the previous and next amino acids in peptides.
        If True, remove them. Default is True.
    remove_modification : bool, optional
        Whether to include modifications in peptides.
        If True, remove them. Default is True.

    Notes
    -----
    The generated features include:
    - mhcflurry_affinity: Binding affinity score
    - mhcflurry_processing_score: Processing score
    - mhcflurry_presentation_score: Presentation score
    - mhcflurry_presentation_percentile: Presentation percentile
    """

    MIN_PEPTIDE_LENGTH = 8
    MAX_PEPTIDE_LENGTH = 15

    def __init__(
        self,
        peptides: List[str],
        alleles: List[str],
        remove_pre_nxt_aa: bool = False,
        remove_modification: bool = True,
        *args,
        **kwargs,
    ):
        self.peptides = peptides
        self.alleles = alleles
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.remove_modification = remove_modification
        self.predictor = Class1PresentationPredictor.load()
        self.predictions = None
        self._raw_predictions = None
        logger.info(
            f"Initialized MHCflurryFeatureGenerator with {len(peptides)} peptides and alleles: {alleles}"
        )

    @property
    def feature_columns(self) -> List[str]:
        """
        Return the list of generated feature column names.

        Returns
        -------
        List[str]
            List of feature column names:
            - mhcflurry_affinity
            - mhcflurry_processing_score
            - mhcflurry_presentation_score
            - mhcflurry_presentation_percentile
        """
        return [
            "mhcflurry_affinity",
            "mhcflurry_processing_score",
            "mhcflurry_presentation_score",
            "mhcflurry_presentation_percentile",
        ]

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
        Preprocess peptide sequence by removing flanking amino acids and modifications.

        Parameters
        ----------
        peptide : str
            Original peptide sequence.

        Returns
        -------
        str
            Preprocessed peptide sequence.
        """
        if self.remove_pre_nxt_aa:
            peptide = utils.strip_flanking_and_charge(peptide)
        if self.remove_modification:
            peptide = utils.remove_modifications(peptide)
        # U -> C
        peptide = peptide.replace("U", "C")
        return peptide

    def _predict(self) -> pd.DataFrame:
        """
        Run MHCflurry predictions using the `predict` method and cache the result.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the prediction results, with the best prediction
            for each sequence.

        Notes
        -----
        This method:
        1. Preprocesses peptides
        2. Filters peptides by length (8-15 amino acids)
        3. Runs MHCflurry predictions
        4. Merges results with original peptides
        """
        if self.predictions is None:
            logger.info("Running MHCflurry predictions.")
            self.predictions = pd.DataFrame(self.peptides, columns=["Peptide"])
            self.predictions["clean_peptide"] = self.predictions["Peptide"].apply(
                self._preprocess_peptides
            )
            # filter out peptides with length between 8 and 15, which are the only lengths supported by MHCflurry
            peptides_to_predict = self.predictions[
                self.predictions["clean_peptide"].apply(
                    lambda x: (
                        MHCflurryFeatureGenerator.MIN_PEPTIDE_LENGTH
                        <= len(x)
                        <= MHCflurryFeatureGenerator.MAX_PEPTIDE_LENGTH
                    )
                )
            ]
            logger.info(
                f"Predicting MHCflurry scores for {len(peptides_to_predict)} peptides. The missing peptides will be filled with median values."
            )
            mhcflurry_results = self.predictor.predict(
                peptides=peptides_to_predict["clean_peptide"].unique().tolist(),
                alleles=self.alleles,
                verbose=0,
            )
            self._raw_predictions = mhcflurry_results.copy()
            self.predictions = self.predictions.merge(
                mhcflurry_results,
                left_on="clean_peptide",
                right_on="peptide",
                how="left",
            )
            self.predictions.drop(columns=["clean_peptide", "peptide"], inplace=True)
        else:
            logger.info("MHCflurry predictions already exist. Skipping prediction.")
        return self.predictions

    @property
    def raw_predictions(self) -> pd.DataFrame:
        """
        Return the raw predictions DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the raw predictions:
            - peptide: Cleaned peptide sequence
            - allele: MHC allele
            - affinity: Binding affinity
            - processing_score: Processing score
            - presentation_score: Presentation score
            - presentation_percentile: Presentation percentile
        """
        if self._raw_predictions is None:
            self._predict()
        return self._raw_predictions

    def get_raw_predictions(self) -> pd.DataFrame:
        """
        Get the raw prediction results DataFrame from MHCflurry.

        Returns
        -------
        pd.DataFrame
            Raw prediction results DataFrame containing:
            - peptide: Cleaned peptide sequence
            - allele: MHC allele
            - affinity: Binding affinity
            - processing_score: Processing score
            - presentation_score: Presentation score
            - presentation_percentile: Presentation percentile
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
        - affinity: Binding affinity
        - processing_score: Processing score
        - presentation_score: Presentation score
        - presentation_percentile: Presentation percentile
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
        Generate MHCflurry features for the provided peptides and alleles.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the peptides and their predicted MHCflurry features:
            - Peptide: Original peptide sequence
            - mhcflurry_affinity: Binding affinity
            - mhcflurry_processing_score: Processing score
            - mhcflurry_presentation_score: Presentation score
            - mhcflurry_presentation_percentile: Presentation percentile

        Notes
        -----
        This method:
        1. Runs MHCflurry predictions
        2. Renames columns to include 'mhcflurry_' prefix
        3. Fills missing values with median values
        4. Returns the final feature DataFrame
        """
        self._predict()
        features_df = self.predictions.copy()
        features_df.rename(
            columns={
                "affinity": "mhcflurry_affinity",
                "processing_score": "mhcflurry_processing_score",
                "presentation_score": "mhcflurry_presentation_score",
                "presentation_percentile": "mhcflurry_presentation_percentile",
            },
            inplace=True,
        )
        features_df.fillna(
            value={
                "mhcflurry_affinity": features_df["mhcflurry_affinity"].median(),
                "mhcflurry_processing_score": features_df["mhcflurry_processing_score"].median(),
                "mhcflurry_presentation_score": features_df[
                    "mhcflurry_presentation_score"
                ].median(),
                "mhcflurry_presentation_percentile": features_df[
                    "mhcflurry_presentation_percentile"
                ].median(),
            },
            inplace=True,
        )
        logger.info(f"Generated MHCflurry features for {len(features_df)} peptides.")
        features_df = features_df[
            [
                "Peptide",
                "mhcflurry_affinity",
                "mhcflurry_processing_score",
                "mhcflurry_presentation_score",
                "mhcflurry_presentation_percentile",
            ]
        ]
        if features_df.isna().sum().sum() > 0:
            logger.warning("NaN values found in the generated features.")
        return features_df

    def get_best_allele(self) -> pd.DataFrame:
        """
        Get the best allele for each peptide.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the best alleles for the peptides:
            - Peptide: Original peptide sequence
            - mhcflurry_best_allele: Best binding allele

        Notes
        -----
        The best allele is determined by the lowest presentation percentile rank.
        """
        best_allele_df = self.predictions[["Peptide", "best_allele"]]
        best_allele_df.rename(columns={"best_allele": "mhcflurry_best_allele"}, inplace=True)

        logger.info(f"Generated best allele information for {len(best_allele_df)} peptides.")

        return best_allele_df

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

    @classmethod
    def from_config(cls, psms, config, params):
        return cls(
            peptides=list(set(psms.peptides)),
            alleles=config.get("allele", []),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
        )


feature_generator_factory.register_generator("MHCflurry", MHCflurryFeatureGenerator)
