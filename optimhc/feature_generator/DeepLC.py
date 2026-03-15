# feature_generator/DeepLC.py
# TODO: Use koina for prediction

import logging
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
from deeplc import DeepLC

from optimhc import utils
from optimhc.feature_generator.base_feature_generator import BaseFeatureGenerator
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)
logging.getLogger("deeplc.feat_extractor").setLevel(logging.CRITICAL)
logging.getLogger("deeplc.feat_extractor").disabled = True


class DeepLCFeatureGenerator(BaseFeatureGenerator):
    """
    Generate DeepLC-based features for rescoring.

    This generator uses DeepLC to predict retention times and calculates various
    features based on the differences between predicted and observed retention times.

    Parameters
    ----------
    psms : PsmContainer
        PSMs to generate features for.
    calibration_criteria_column : str
        Column name in the PSMs DataFrame to use for DeepLC calibration.
    lower_score_is_better : bool, optional
        Whether a lower PSM score denotes a better matching PSM. Default is False.
    calibration_set_size : int or float, optional
        Amount of best PSMs to use for DeepLC calibration. If this value is lower
        than the number of available PSMs, all PSMs will be used. Default is 0.15.
    processes : int, optional
        Number of processes to use in DeepLC. Default is 1.
    model_path : str, optional
        Path to the DeepLC model. If None, the default model will be used.
    remove_pre_nxt_aa : bool, optional
        Whether to remove the first and last amino acids from the peptide sequence.
        Default is True.
    mod_dict : dict, optional
        Dictionary of modifications to be used for DeepLC. If None, no modifications
        will be used.

    Notes
    -----
    DeepLC retraining is on by default. Add ``deeplc_retrain: False`` as a keyword
    argument to disable retraining.

    The generated features include:
    - observed_retention_time: Original retention time from the data
    - predicted_retention_time: DeepLC predicted retention time
    - retention_time_diff: Difference between predicted and observed times
    - abs_retention_time_diff: Absolute difference between predicted and observed times
    - retention_time_ratio: Ratio of min(pred,obs) to max(pred,obs)
    """

    def __init__(
        self,
        psms: PsmContainer,
        calibration_criteria_column: str,
        lower_score_is_better: bool = False,
        calibration_set_size: Union[int, float, None] = None,
        processes: int = 1,
        model_path: Optional[str] = None,
        remove_pre_nxt_aa: bool = True,
        mod_dict: Optional[Dict[str, str]] = None,
        *args,
        **kwargs,
    ):
        """
        Generate DeepLC-based features for rescoring.

        DeepLC retraining is on by default. Add ``deeplc_retrain: False`` as a keyword argument to
        disable retraining.

        Parameters
        ----------
        psms : PsmContainer
            PSMs to generate features for.
        calibration_criteria_column : str
            Column name in the PSMs DataFrame to use for DeepLC calibration.
        lower_score_is_better : bool
            Whether a lower PSM score denotes a better matching PSM. Default: False.
        calibration_set_size : int or float
            Amount of best PSMs to use for DeepLC calibration. If this value is lower
            than the number of available PSMs, all PSMs will be used. (default: 0.15)
        processes : int or None
            Number of processes to use in DeepLC. Defaults to 1.
        model_path : str
            Path to the DeepLC model. If None, the default model will be used.
        remove_pre_nxt_aa : bool
            Whether to remove the first and last amino acids from the peptide sequence.
            Default: True.
        mod_dict : dict
            Dictionary of modifications to be used for DeepLC. If None, no modifications
            will be used.
        *args : list
            Additional positional arguments are passed to DeepLC.
        **kwargs : dict
            Additional keyword arguments are passed to DeepLC.
        """
        self.psms = psms
        self.lower_score_is_better = lower_score_is_better
        self.calibration_criteria_column = calibration_criteria_column
        self.calibration_set_size = calibration_set_size
        self.processes = processes
        self.model_path = model_path
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.mod_dict = mod_dict
        self.deeplc_df = self._get_deeplc_df()
        self.DeepLC = DeepLC
        self._raw_predictions = None
        if model_path is not None:
            self.deeplc_predictor = self.DeepLC(
                n_jobs=self.processes,
                path_model=model_path,
            )
        else:
            self.deeplc_predictor = self.DeepLC(n_jobs=self.processes)
        logger.info(
            f"Initialized DeepLCFeatureGenerator with {len(self.psms)} PSMs."
            f" Calibration criteria: {self.calibration_criteria_column}."
            f" Lower score is better: {self.lower_score_is_better}."
            f" Calibration set size: {self.calibration_set_size}."
            f" Processes: {self.processes}."
            f" Model path: {self.model_path}."
        )

    @property
    def feature_columns(self) -> List[str]:
        """
        Return the list of generated feature column names.

        Returns
        -------
        List[str]
            List of feature column names:
            - observed_retention_time
            - predicted_retention_time
            - retention_time_diff
            - abs_retention_time_diff
            - retention_time_ratio
        """
        return [
            "observed_retention_time",
            "predicted_retention_time",
            "retention_time_diff",
            "abs_retention_time_diff",
            "retention_time_ratio",
        ]

    @property
    def id_column(self) -> List[str]:
        """
        Return the list of input columns required for the feature generator.

        Returns
        -------
        List[str]
            List of input columns required for feature generation.
            Currently returns an empty list as the required columns are
            handled internally by the PsmContainer.
        """
        return [""]

    def _get_deeplc_df(self):
        """
        Extract the format required by DeepLC, while retaining necessary original information.

        Returns
        -------
        pd.DataFrame
            DataFrame with the required DeepLC format and original information:
            - original_seq: Original peptide sequence
            - label: Target/decoy label
            - seq: Cleaned peptide sequence
            - modifications: Unimod format modifications
            - tr: Retention time
            - score: Calibration criteria score

        Raises
        ------
        ValueError
            If retention time column is not found in the PSMs DataFrame.

        Notes
        -----
        This method prepares the data in the format required by DeepLC,
        including cleaning peptide sequences and converting modifications
        to Unimod format.
        """
        df_deeplc = pd.DataFrame()
        df_psm = self.psms.psms
        df_deeplc["original_seq"] = df_psm[self.psms.peptide_column]
        df_deeplc["label"] = df_psm[self.psms.label_column]

        if self.remove_pre_nxt_aa:
            df_deeplc["seq"] = df_deeplc["original_seq"].apply(utils.strip_flanking_and_charge)
        else:
            df_deeplc["seq"] = df_deeplc["original_seq"]

        # Apply extract_unimod_from_peptidoform once and store both results.
        if self.mod_dict is None:
            logger.warning("No mod_dict provided. Removing modifications.")
            df_deeplc["seq"] = df_deeplc["seq"].apply(lambda x: utils.remove_modifications(x))
            df_deeplc["modifications"] = ""
        else:
            extracted_results = df_deeplc["seq"].apply(
                lambda x: utils.extract_unimod_from_peptidoform(x, mod_dict=self.mod_dict)
            )
            df_deeplc["seq"] = extracted_results.apply(lambda x: x[0])
            df_deeplc["modifications"] = extracted_results.apply(lambda x: x[1])

        if self.psms.retention_time_column is None:
            raise ValueError("DeepLC requires retention time values.")

        df_deeplc["tr"] = df_psm[self.psms.retention_time_column]
        df_deeplc["score"] = df_psm[self.calibration_criteria_column]

        logger.debug("DeepLC input DataFrame:")
        logger.debug(df_deeplc)

        return df_deeplc

    def generate_features(self) -> pd.DataFrame:
        """
        Generate DeepLC features for the provided PSMs.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the PSMs with added DeepLC features:
            - original_seq: Original peptide sequence
            - observed_retention_time: Original retention time
            - predicted_retention_time: DeepLC predicted retention time
            - retention_time_diff: Difference between predicted and observed times
            - abs_retention_time_diff: Absolute difference between predicted and observed times
            - retention_time_ratio: Ratio of min(pred,obs) to max(pred,obs)

        Notes
        -----
        This method:
        1. Prepares data in DeepLC format
        2. Calibrates DeepLC if calibration set is specified
        3. Predicts retention times
        4. Calculates various retention time-based features
        5. Handles missing values by imputing with median values
        """
        logger.info("Generating DeepLC features.")

        # Extract DeepLC input DataFrame
        self.deeplc_df = self._get_deeplc_df()

        # Calibrate DeepLC predictor
        if self.calibration_set_size:
            calibration_df = self._get_calibration_psms(self.deeplc_df)
            logger.debug(f"Calibrating DeepLC with {len(calibration_df)} PSMs.")
            self.deeplc_predictor.calibrate_preds(
                seq_df=calibration_df[["seq", "tr", "modifications"]]
            )

        # Predict retention times
        logger.info("Predicting retention times using DeepLC.")
        predictions = self.deeplc_predictor.make_preds(
            seq_df=self.deeplc_df[["seq", "tr", "modifications"]]
        )

        self._raw_predictions = pd.DataFrame(
            {
                "peptide": self.deeplc_df["seq"],
                "predicted_rt": predictions,
                "observed_rt": self.deeplc_df["tr"],
                "modifications": self.deeplc_df["modifications"],
            }
        )

        # Calculate retention time differences
        rt_diffs = predictions - self.deeplc_df["tr"]
        self.deeplc_df["predicted_retention_time"] = predictions
        self.deeplc_df["retention_time_diff"] = rt_diffs

        result_df = pd.DataFrame()
        result_df["original_seq"] = self.deeplc_df["original_seq"]
        result_df["observed_retention_time"] = self.deeplc_df["tr"]
        result_df["predicted_retention_time"] = self.deeplc_df["predicted_retention_time"]
        result_df["retention_time_diff"] = self.deeplc_df["retention_time_diff"]
        result_df["abs_retention_time_diff"] = self.deeplc_df["retention_time_diff"].abs()

        # Adopt from 'DeepRescore2': RTR = min(pred, obs) / max(pred, obs)
        result_df["retention_time_ratio"] = np.minimum(
            result_df["predicted_retention_time"], result_df["observed_retention_time"]
        ) / np.maximum(result_df["predicted_retention_time"], result_df["observed_retention_time"])

        for col in self.feature_columns:
            nan_rows = result_df[result_df[col].isna()]
            if not nan_rows.empty:
                logger.warning(
                    f"Column {col} contains NaN values. Rows with NaN values:\n{nan_rows}"
                )
            median_value = result_df[col].median()
            result_df[col].fillna(median_value, inplace=True)
            result_df[col] = result_df[col].astype(float)

        return result_df

    def _get_calibration_psms(self, deeplc_df: pd.DataFrame) -> pd.DataFrame:
        """
        Get the best scoring PSMs for calibration based on the calibration criteria.

        Parameters
        ----------
        deeplc_df : pd.DataFrame
            DataFrame containing DeepLC input data.

        Returns
        -------
        pd.DataFrame
            DataFrame of PSMs selected for calibration, containing only target PSMs.

        Raises
        ------
        ValueError
            If calibration_set_size is a float not between 0 and 1.
        TypeError
            If calibration_set_size is neither int nor float.

        Notes
        -----
        This method:
        1. Sorts PSMs based on calibration criteria
        2. Selects top N PSMs based on calibration_set_size
        3. Filters to keep only target PSMs
        """
        logger.debug("Selecting PSMs for calibration.")

        # Sort PSMs based on calibration criteria
        sorted_psms = deeplc_df.sort_values(by="score", ascending=self.lower_score_is_better)

        # Select calibration set
        if isinstance(self.calibration_set_size, float):
            if not 0 < self.calibration_set_size <= 1:
                logger.error("calibration_set_size as float must be between 0 and 1.")
                raise ValueError(
                    "If `calibration_set_size` is a float, it must be between 0 and 1."
                )
            n_cal = int(len(sorted_psms) * self.calibration_set_size)
        elif isinstance(self.calibration_set_size, int):
            n_cal = self.calibration_set_size
            if n_cal > len(sorted_psms):
                logger.warning(
                    f"Requested calibration_set_size ({n_cal}) exceeds number of PSMs ({len(sorted_psms)}). Using all PSMs for calibration."
                )
                n_cal = len(sorted_psms)
        else:
            logger.error("calibration_set_size must be either int or float.")
            raise TypeError(
                f"Expected int or float for `calibration_set_size`. Got {type(self.calibration_set_size)} instead."
            )

        calibration_psms = sorted_psms.head(n_cal)
        logger.debug(f"Selected {n_cal} PSMs for calibration.")
        calibration_psms = calibration_psms[calibration_psms["label"]]
        logger.debug(f"Selected {len(calibration_psms)} target PSMs for calibration.")
        return calibration_psms

    def get_full_data(self) -> pd.DataFrame:
        """
        Get the full DeepLC DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the DeepLC input data with all columns:
            - original_seq: Original peptide sequence
            - label: Target/decoy label
            - seq: Cleaned peptide sequence
            - modifications: Unimod format modifications
            - tr: Retention time
            - score: Calibration criteria score
            - predicted_retention_time: DeepLC predicted retention time
            - retention_time_diff: Difference between predicted and observed times
        """
        return self.deeplc_df

    @property
    def raw_predictions(self) -> pd.DataFrame:
        """
        Get the raw predictions DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the raw predictions:
            - peptide: Cleaned peptide sequence
            - predicted_rt: DeepLC predicted retention time
            - observed_rt: Original retention time
            - modifications: Unimod format modifications

        Notes
        -----
        If predictions haven't been generated yet, this will trigger
        feature generation automatically.
        """
        if self._raw_predictions is None:
            self.generate_features()
        return self._raw_predictions

    def get_raw_predictions(self) -> pd.DataFrame:
        """
        Get the raw predictions DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing the raw predictions:
            - peptide: Cleaned peptide sequence
            - predicted_rt: DeepLC predicted retention time
            - observed_rt: Original retention time
            - modifications: Unimod format modifications

        Notes
        -----
        This is a convenience method that returns the same data as the
        raw_predictions property.
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
        - predicted_rt: DeepLC predicted retention time
        - observed_rt: Original retention time
        - modifications: Unimod format modifications
        """
        if "index" not in kwargs:
            kwargs["index"] = False
        if self.raw_predictions is not None:
            self.raw_predictions.to_csv(file_path, **kwargs)
            logger.info(f"Raw predictions saved to {file_path}")
        else:
            logger.warning("Raw predictions have not been generated yet.")
