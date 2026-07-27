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

    MHCflurry's predict() API returns per-peptide best-allele results directly,
    so this generator always operates in "best" mode with no per-allele breakdown.

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC allele names (e.g., ['HLA-A*02:02', 'HLA-B*07:02']).
    remove_pre_nxt_aa : bool, optional
        Whether to remove flanking previous/next amino acids from peptide sequences.
        Default is False.
    remove_modification : bool, optional
        Whether to remove peptide modification annotations. Default is True.

    Notes
    -----
    Generated features:
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
        return [
            "mhcflurry_affinity",
            "mhcflurry_processing_score",
            "mhcflurry_presentation_score",
            "mhcflurry_presentation_percentile",
        ]

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

    def _predict(self) -> pd.DataFrame:
        if self.predictions is not None:
            logger.info("MHCflurry predictions already exist. Skipping prediction.")
            return self.predictions

        logger.info("Running MHCflurry predictions.")
        # Use a local variable so self.predictions is only set on full success,
        # preventing a partial state from poisoning the cache if an exception occurs.
        predictions = pd.DataFrame(self.peptides, columns=["Peptide"])
        predictions["clean_peptide"] = predictions["Peptide"].apply(self._preprocess_peptides)

        peptides_to_predict = predictions[
            predictions["clean_peptide"].apply(
                lambda x: (
                    MHCflurryFeatureGenerator.MIN_PEPTIDE_LENGTH
                    <= len(x)
                    <= MHCflurryFeatureGenerator.MAX_PEPTIDE_LENGTH
                )
            )
        ]
        logger.info(
            f"Predicting MHCflurry scores for {len(peptides_to_predict)} peptides. "
            "Missing peptides will be filled with median values."
        )

        mhcflurry_results = self.predictor.predict(
            peptides=peptides_to_predict["clean_peptide"].unique().tolist(),
            alleles=self.alleles,
            verbose=0,
        )
        self._raw_predictions = mhcflurry_results.copy()

        predictions = predictions.merge(
            mhcflurry_results,
            left_on="clean_peptide",
            right_on="peptide",
            how="left",
        )
        predictions.drop(columns=["clean_peptide", "peptide"], inplace=True)

        self.predictions = predictions
        return self.predictions

    @property
    def raw_predictions(self) -> pd.DataFrame:
        if self._raw_predictions is None:
            self._predict()
        return self._raw_predictions

    def generate_features(self) -> pd.DataFrame:
        """
        Generate MHCflurry features for the provided peptides and alleles.

        Returns a DataFrame with per-peptide best-allele predictions.
        Missing values are filled with column medians.
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

        for col in self.feature_columns:
            if col in features_df.columns:
                features_df[col] = features_df[col].fillna(features_df[col].median())

        features_df = features_df[["Peptide"] + self.feature_columns]

        if features_df.isna().sum().sum() > 0:
            logger.warning("NaN values found in the generated features.")
        logger.info(f"Generated MHCflurry features for {len(features_df)} peptides.")
        return features_df

    @classmethod
    def from_config(cls, psms, config, params):
        return cls(
            peptides=list(set(psms.df["sequence"])),
            alleles=config.get("allele", []),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
        )


feature_generator_factory.register_generator("MHCflurry", MHCflurryFeatureGenerator)
