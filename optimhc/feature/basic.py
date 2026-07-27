import logging
from typing import List

import pandas as pd
from scipy.stats import entropy  # Import entropy from scipy

from optimhc import utils
from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)


class BasicFeatureGenerator(BaseFeatureGenerator):
    """
    Feature generator that generates basic features from peptide sequences.

    This generator calculates features such as peptide length, proportion of unique amino acids,
    Shannon entropy of amino acid distribution, difference between peptide length and average peptide length,
    and count of unique amino acids.

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences to generate features for.
    remove_pre_nxt_aa : bool, optional
        Whether to remove flanking previous/next amino acids from peptide sequences.
        Default is True.
    remove_modification : bool, optional
        Whether to remove peptide modification annotations. Default is True.

    Notes
    -----
    The generated features include:
    - length_diff_from_avg: Difference between peptide length and average length
    - abs_length_diff_from_avg: Absolute difference between peptide length and average length
    - unique_aa_count: Number of unique amino acids in the peptide
    - unique_aa_proportion: Proportion of unique amino acids in the peptide
    - shannon_entropy: Shannon entropy of amino acid distribution
    """

    def __init__(
        self,
        peptides: List[str],
        remove_pre_nxt_aa: bool = True,
        remove_modification: bool = True,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.peptides = peptides
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.remove_modification = remove_modification
        self.avg_length = None
        logger.info(f"Initialized BasicFeatureGenerator with {len(peptides)} peptides.")

    @property
    def feature_columns(self) -> List[str]:
        """Return the list of generated feature column names."""
        return [
            #'peptide_length',
            "length_diff_from_avg",
            "abs_length_diff_from_avg",
            "unique_aa_count",
            "unique_aa_proportion",
            "shannon_entropy",
        ]

    @property
    def id_column(self) -> List[str]:
        """
        Return the list of input columns required for feature generation.

        Returns
        -------
        List[str]
            List of input column names required for feature generation.
            Currently only requires 'Peptide' column.
        """
        return ["Peptide"]

    def _preprocess_peptide(self, peptide: str) -> str:
        """
        Preprocess peptide sequence by removing adjacent amino acids and modifications.

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
        return peptide

    def _shannon_entropy(self, sequence: str) -> float:
        """
        Calculate the Shannon entropy of a peptide sequence.

        Parameters
        ----------
        sequence : str
            Peptide sequence.

        Returns
        -------
        float
            Shannon entropy value.
        """
        if len(sequence) == 0:
            return 0.0
        # Calculate frequency of each unique amino acid
        bases = list(set(sequence))
        freq_list = [sequence.count(base) / len(sequence) for base in bases]
        return entropy(freq_list, base=2)

    def generate_features(self) -> pd.DataFrame:
        """
        Generate basic features for the provided peptides.

        Returns
        -------
        pd.DataFrame
            DataFrame containing peptides and their computed features:
            - length_diff_from_avg: Difference from average peptide length
            - abs_length_diff_from_avg: Absolute difference from average length
            - unique_aa_count: Number of unique amino acids
            - unique_aa_proportion: Proportion of unique amino acids
            - shannon_entropy: Shannon entropy of amino acid distribution

        Raises
        ------
        ValueError
            If NaN values are found in the generated features.

        Notes
        -----
        All features are converted to float type before returning.
        The method calculates average peptide length across all peptides
        and uses it as a reference for length-based features.
        """
        logger.info("Generating basic features.")
        peptides_df = pd.DataFrame(self.peptides, columns=["Peptide"])
        peptides_df["clean_peptide"] = peptides_df["Peptide"].apply(self._preprocess_peptide)
        peptides_df["peptide_length"] = peptides_df["clean_peptide"].apply(len)
        self.avg_length = peptides_df["peptide_length"].mean()
        peptides_df["length_diff_from_avg"] = peptides_df["peptide_length"] - self.avg_length
        peptides_df["abs_length_diff_from_avg"] = peptides_df["length_diff_from_avg"].abs()
        peptides_df["unique_aa_count"] = peptides_df["clean_peptide"].apply(lambda x: len(set(x)))
        peptides_df["unique_aa_proportion"] = (
            peptides_df["unique_aa_count"] / peptides_df["peptide_length"]
        )
        peptides_df["shannon_entropy"] = peptides_df["clean_peptide"].apply(self._shannon_entropy)
        features_df = peptides_df[["Peptide"] + self.feature_columns]
        # Fix SettingWithCopyWarning: make an explicit copy before assignment
        features_df = features_df.copy()
        for col in self.feature_columns:
            features_df[col] = features_df[col].astype(float)
        if features_df.isna().sum().sum() > 0:
            logger.error("NaN values found in the generated features.")
            raise ValueError("NaN values found in the generated features.")

        logger.info(f"Generated basic features for {len(features_df)} peptides.")
        return features_df

    @classmethod
    def from_config(cls, psms, config, params):
        return cls(
            peptides=psms.df["sequence"].tolist(),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
        )

    def apply(self, psms):
        features = self.generate_features()
        features.insert(0, "psm_id", psms.df["psm_id"].to_numpy())
        psms.add_features(features, on="psm_id", columns=self.feature_columns)


feature_generator_factory.register_generator("Basic", BasicFeatureGenerator)
