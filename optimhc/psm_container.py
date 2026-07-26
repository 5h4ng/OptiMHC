"""OptiMHC peptide-spectrum match table."""

from __future__ import annotations

import numpy as np
import pandas as pd

REQUIRED_COLUMNS = (
    "psm_id",
    "run",
    "scan",
    "rank",
    "sequence",
    "mods",
    "mod_sites",
    "charge",
    "proteins",
    "is_decoy",
)


class PsmContainer:
    """Store PSM rows and the columns used for rescoring.

    The DataFrame is stored by reference and must contain ``REQUIRED_COLUMNS``.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        feature_columns: list[str] | tuple[str, ...] = (),
    ) -> None:
        self.df = df
        feature_names = tuple(feature_columns)
        self.feature_columns = (
            feature_names if "rank" in feature_names else (*feature_names, "rank")
        )
        self._validate()

    def _validate(self) -> None:
        missing = [column for column in REQUIRED_COLUMNS if column not in self.df.columns]
        if missing:
            raise ValueError(f"Missing required PSM columns: {missing}")

        missing_features = [
            column for column in self.feature_columns if column not in self.df.columns
        ]
        if missing_features:
            raise ValueError(f"Feature columns not found in PSM data: {missing_features}")
        if len(set(self.feature_columns)) != len(self.feature_columns):
            raise ValueError("Feature columns must be unique.")
        if self.df.loc[:, REQUIRED_COLUMNS].isna().any().any():
            raise ValueError("Required PSM columns cannot contain missing values.")
        if not pd.api.types.is_integer_dtype(self.df["psm_id"]):
            raise ValueError("Column 'psm_id' must contain integers.")
        if not self.df["psm_id"].is_unique:
            raise ValueError("Column 'psm_id' must be unique.")
        if not pd.api.types.is_bool_dtype(self.df["is_decoy"]):
            raise ValueError("Column 'is_decoy' must be boolean.")

        mismatched_mods = [
            index
            for index, (mods, sites) in enumerate(zip(self.df["mods"], self.df["mod_sites"]))
            if _token_count(mods) != _token_count(sites)
        ]
        if mismatched_mods:
            raise ValueError("Columns 'mods' and 'mod_sites' must contain aligned values.")

    def __len__(self) -> int:
        return len(self.df)

    def __repr__(self) -> str:
        n_target = int((~self.df["is_decoy"]).sum())
        n_decoy = int(self.df["is_decoy"].sum())
        n_spectra = self.df[["run", "scan"]].drop_duplicates().shape[0]
        n_runs = self.df["run"].nunique()
        return (
            f"PsmContainer({len(self)} PSMs ({n_target} target, {n_decoy} decoy), "
            f"{n_spectra} spectra, {n_runs} run(s), "
            f"{len(self.feature_columns)} features)"
        )

    def add_features(
        self,
        features: pd.DataFrame,
        *,
        on: str | list[str] | tuple[str, ...],
        columns: list[str] | tuple[str, ...],
    ) -> None:
        """Add numeric feature columns by matching the specified key columns.

        The feature keys must be unique and cover every PSM key. New column
        names must not already exist, and all feature values must be finite
        numbers.
        """
        keys = [on] if isinstance(on, str) else list(on)
        new_columns = list(columns)

        if not keys or not new_columns:
            raise ValueError("Feature join keys and columns must be explicit.")
        missing_left = [column for column in keys if column not in self.df.columns]
        missing_right = [
            column for column in (*keys, *new_columns) if column not in features.columns
        ]
        if missing_left or missing_right:
            raise ValueError(
                f"Missing feature join columns: PSM={missing_left}, features={missing_right}"
            )
        conflicts = sorted(set(new_columns).intersection(self.df.columns))
        if conflicts:
            raise ValueError(f"Feature columns already exist in PSM data: {conflicts}")
        if features.duplicated(keys).any():
            raise ValueError(f"Feature keys must be unique: {keys}")

        key_coverage = (
            self.df[keys]
            .drop_duplicates()
            .merge(features[keys].drop_duplicates(), how="outer", on=keys, indicator=True)
        )
        if not key_coverage["_merge"].eq("both").all():
            raise ValueError("Feature keys must exactly cover the PSM keys.")

        numeric_values = features[new_columns].apply(pd.to_numeric, errors="raise")
        if not np.isfinite(numeric_values.to_numpy(dtype=float)).all():
            raise ValueError("Feature columns must contain only finite numeric values.")

        feature_df = features[[*keys, *new_columns]].copy()
        feature_df[new_columns] = numeric_values
        matched_features = self.df[keys].merge(
            feature_df,
            how="left",
            on=keys,
            sort=False,
            validate="many_to_one",
        )
        for column in new_columns:
            self.df[column] = matched_features[column].to_numpy()
        self.feature_columns = (*self.feature_columns, *new_columns)


def _token_count(value: object) -> int:
    text = str(value)
    return 0 if not text else len(text.split(";"))
