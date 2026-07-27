import logging
import os
from multiprocessing import Pool, cpu_count
from typing import List, Optional

import pandas as pd
from mhctools import NetMHCpan41_BA
from tqdm import tqdm

from optimhc import utils
from optimhc.feature.factory import feature_generator_factory

from .base_feature_generator import BaseFeatureGenerator

logger = logging.getLogger(__name__)


# Each worker process gets its own copy of this global.
# Pool(initializer=_init_worker) sets it once per process,
# so netMHCpan -listMHC is called O(processes) not O(chunks).
_worker_predictor: Optional[NetMHCpan41_BA] = None


def _init_worker(alleles: List[str], program_name: str) -> None:
    global _worker_predictor
    # mhctools reuses program_name as a temporary-file prefix, so it must not
    # contain directory separators.
    #
    # Keep the executable discoverable by adding its directory to this process's
    # PATH, then pass only the basename to mhctools.
    exe_dir = os.path.dirname(program_name)
    if exe_dir:
        os.environ["PATH"] = exe_dir + os.pathsep + os.environ.get("PATH", "")
        program_name = os.path.basename(program_name)
    _worker_predictor = NetMHCpan41_BA(alleles=alleles, program_name=program_name)


def _predict_peptide_chunk(peptides_chunk: List[str]) -> pd.DataFrame:
    return _worker_predictor.predict_peptides(peptides_chunk).to_dataframe()


class NetMHCpanFeatureGenerator(BaseFeatureGenerator):
    """
    Generate NetMHCpan binding prediction features for MHC class I alleles.

    For each peptide, reports the prediction for the allele with the lowest
    percentile rank (strongest predicted binder).

    Parameters
    ----------
    peptides : List[str]
        List of peptide sequences.
    alleles : List[str]
        List of MHC allele names (e.g., ['HLA-A*02:01', 'HLA-B*07:02']).
    remove_pre_nxt_aa : bool, optional
        Whether to remove flanking previous/next amino acids from peptide sequences.
        Default is False.
    remove_modification : bool, optional
        Whether to remove peptide modification annotations. Default is True.
    n_processes : int, optional
        Number of processes to use for multiprocessing. Default is 1.
    show_progress : bool, optional
        Whether to display a progress bar. Default is False.
    executable_path : str, optional
        Path to the netMHCpan executable. Defaults to "netMHCpan" (PATH lookup).

    Notes
    -----
    Generated features:
    - netmhcpan_score: Raw binding score
    - netmhcpan_affinity: Binding affinity in nM
    - netmhcpan_percentile_rank: Percentile rank
    """

    MIN_PEPTIDE_LENGTH = 8
    MAX_PEPTIDE_LENGTH = 30
    CHUNKSIZE = 250

    def __init__(
        self,
        peptides: List[str],
        alleles: List[str],
        remove_pre_nxt_aa: bool = False,
        remove_modification: bool = True,
        n_processes: int = 1,
        show_progress: bool = False,
        executable_path: str = "netMHCpan",
    ):
        self.peptides = peptides
        self.alleles = alleles
        self.remove_pre_nxt_aa = remove_pre_nxt_aa
        self.remove_modification = remove_modification
        self.n_processes = min(n_processes, cpu_count())
        self.show_progress = show_progress
        self.executable_path = executable_path

        # mhctools reuses program_name as a temporary-file prefix, so it must not
        # contain directory separators.
        #
        # Keep the executable discoverable by adding its directory to this process's
        # PATH, then pass only the basename to mhctools.
        exe_dir = os.path.dirname(executable_path)
        if exe_dir:
            os.environ["PATH"] = exe_dir + os.pathsep + os.environ.get("PATH", "")
            executable_path = os.path.basename(executable_path)

        # Only create a predictor when running single-process; multiprocessing uses
        # _init_worker to create one per worker process instead.
        self.predictor = (
            NetMHCpan41_BA(alleles=self.alleles, program_name=executable_path)
            if self.n_processes == 1
            else None
        )
        self.predictions = None
        self._raw_predictions = None
        logger.info(
            f"Initialized NetMHCpanFeatureGenerator with {len(peptides)} peptides, "
            f"alleles: {alleles}, n_processes: {self.n_processes}"
        )

    @property
    def feature_columns(self) -> List[str]:
        return [
            "netmhcpan_score",
            "netmhcpan_affinity",
            "netmhcpan_percentile_rank",
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

    def _predict_multiprocessing(self, peptides_to_predict: List[str]) -> pd.DataFrame:
        """
        Run NetMHCpan predictions using multiprocessing.

        One predictor is initialized per worker process (via Pool initializer),
        so allele metadata is loaded once per process rather than once per chunk.
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
        # Use a local variable so self.predictions is only set on full success,
        # preventing a partial state from poisoning the cache if an exception occurs.
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

        For each peptide, returns the prediction for the allele with the lowest
        percentile rank. Missing values are filled with column medians.
        """
        predictions_df = self._predict()
        features_df = pd.DataFrame({"Peptide": self.peptides})
        features_df = self._select_allele_features(predictions_df, features_df)
        features_df = self._fill_missing_values(features_df)

        features_df = features_df[["Peptide"] + self.feature_columns]

        if features_df.isna().sum().sum() > 0:
            logger.warning(
                "NaN values still exist in the generated features after filling with median."
            )
        return features_df

    def _select_allele_features(
        self, predictions_df: pd.DataFrame, features_df: pd.DataFrame
    ) -> pd.DataFrame:
        logger.info("Selecting per-peptide allele predictions.")
        valid_predictions = predictions_df.dropna(subset=["percentile_rank"])

        if valid_predictions.empty:
            logger.warning("No valid predictions available.")
            allele_features = pd.DataFrame(
                {
                    "Peptide": self.peptides,
                    "netmhcpan_allele": ["Unknown"] * len(self.peptides),
                    "netmhcpan_score": [pd.NA] * len(self.peptides),
                    "netmhcpan_affinity": [pd.NA] * len(self.peptides),
                    "netmhcpan_percentile_rank": [pd.NA] * len(self.peptides),
                }
            )
        else:
            idx = valid_predictions.groupby("Peptide")["percentile_rank"].idxmin()
            allele_features = valid_predictions.loc[idx].rename(
                columns={
                    "allele": "netmhcpan_allele",
                    "score": "netmhcpan_score",
                    "affinity": "netmhcpan_affinity",
                    "percentile_rank": "netmhcpan_percentile_rank",
                }
            )[
                [
                    "Peptide",
                    "netmhcpan_allele",
                    "netmhcpan_score",
                    "netmhcpan_affinity",
                    "netmhcpan_percentile_rank",
                ]
            ]

            missing_peptides = set(self.peptides) - set(allele_features["Peptide"])
            if missing_peptides:
                logger.warning(f"Found {len(missing_peptides)} peptides with no prediction.")
                missing_features = pd.DataFrame(
                    {
                        "Peptide": list(missing_peptides),
                        "netmhcpan_allele": ["Unknown"] * len(missing_peptides),
                        "netmhcpan_score": [pd.NA] * len(missing_peptides),
                        "netmhcpan_affinity": [pd.NA] * len(missing_peptides),
                        "netmhcpan_percentile_rank": [pd.NA] * len(missing_peptides),
                    }
                )
                allele_features = pd.concat([allele_features, missing_features], ignore_index=True)

        features_df = features_df.merge(allele_features, on="Peptide", how="left")
        return features_df

    def _fill_missing_values(self, features_df: pd.DataFrame) -> pd.DataFrame:
        if "netmhcpan_allele" in features_df.columns:
            features_df["netmhcpan_allele"] = features_df["netmhcpan_allele"].fillna("Unknown")

        for col in ["netmhcpan_score", "netmhcpan_affinity", "netmhcpan_percentile_rank"]:
            if col in features_df.columns and features_df[col].isna().any():
                median_value = features_df[col].median()
                features_df[col] = features_df[col].fillna(
                    0 if pd.isna(median_value) else median_value
                )
        return features_df

    @classmethod
    def from_config(cls, psms, config, params):
        return cls(
            peptides=list(set(psms.df["sequence"])),
            alleles=config.get("allele", []),
            remove_pre_nxt_aa=config["removePreNxtAA"],
            remove_modification=True,
            n_processes=config.get("numProcesses", 1),
            show_progress=config.get("showProgress", False),
            executable_path=params.get("executablePath", "netMHCpan"),
        )


feature_generator_factory.register_generator("NetMHCpan", NetMHCpanFeatureGenerator)
