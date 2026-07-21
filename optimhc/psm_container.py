"""Canonical peptide-spectrum match table."""

from __future__ import annotations

from collections.abc import Iterable, Sequence

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
    """A small mutable wrapper around the canonical PSM DataFrame."""

    def __init__(
        self,
        df: pd.DataFrame,
        feature_columns: Iterable[str] = (),
    ) -> None:
        self.df = df
        declared = tuple(feature_columns)
        self.feature_columns = declared if "rank" in declared else (*declared, "rank")
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
        return (
            f"PsmContainer({len(self)} candidates, "
            f"{self.df[['run', 'scan']].drop_duplicates().shape[0]} spectra, "
            f"{len(self.feature_columns)} features)"
        )

    def add_features(
        self,
        features: pd.DataFrame,
        *,
        on: str | Sequence[str],
        columns: Sequence[str],
    ) -> None:
        """Attach explicitly declared numeric features using canonical keys."""
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

        key_coverage = self.df[keys].drop_duplicates().merge(
            features[keys].drop_duplicates(), how="outer", on=keys, indicator=True
        )
        if not key_coverage["_merge"].eq("both").all():
            raise ValueError("Feature keys must exactly cover the PSM keys.")

        numeric = features[new_columns].apply(pd.to_numeric, errors="raise")
        if not np.isfinite(numeric.to_numpy(dtype=float)).all():
            raise ValueError("Feature columns must contain only finite numeric values.")

        normalized = features[[*keys, *new_columns]].copy()
        normalized[new_columns] = numeric
        merged = self.df[keys].merge(
            normalized,
            how="left",
            on=keys,
            sort=False,
            validate="many_to_one",
        )
        for column in new_columns:
            self.df[column] = merged[column].to_numpy()
        self.feature_columns = (*self.feature_columns, *new_columns)


def _token_count(value: object) -> int:
    text = str(value)
    return 0 if not text else len(text.split(";"))
