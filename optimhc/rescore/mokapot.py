"""Convert OptiMHC PSMs to Mokapot input and run rescoring."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import mokapot as mokapot_lib
import numpy as np
import pandas as pd

from optimhc.peptidoform import to_proforma
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def mokapot_spec_id(run: object, scan: object, charge: object) -> str:
    """Build a Mokapot SpecId from run, scan, and charge.

    Candidates with the same run, scan, and charge intentionally share a
    SpecId. Mokapot groups spectra by ``filename`` and ``ScanNr``.
    """
    return f"{run}.{int(scan)}.{int(scan)}.{int(charge)}"


def to_mokapot_dataframe(
    psms: PsmContainer,
    feature_columns: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Convert PSMs to a PIN-format DataFrame.

    Requested features must be declared by ``psms`` and contain only finite
    numeric values. ``rank`` is always included as a feature. If no charge
    feature is selected, the integer ``charge`` column is added as ``Charge``.
    """
    selected_columns = _selected_feature_columns(psms, feature_columns)
    unknown_columns = [column for column in selected_columns if column not in psms.feature_columns]
    if unknown_columns:
        raise ValueError(f"Unknown rescoring features: {unknown_columns}")

    psm_df = psms.df
    try:
        feature_values = psm_df.loc[:, selected_columns].apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise ValueError("Rescoring features must contain only finite numeric values.") from error
    if not np.isfinite(feature_values.to_numpy(dtype=float)).all():
        raise ValueError("Rescoring features must contain only finite numeric values.")
    pin_df = pd.DataFrame(
        {
            "SpecId": [
                mokapot_spec_id(run, scan, charge)
                for run, scan, charge in zip(psm_df["run"], psm_df["scan"], psm_df["charge"])
            ],
            "Label": psm_df["is_decoy"].map({False: 1, True: -1}).astype(int),
            "ScanNr": psm_df["scan"].astype(int),
            "filename": psm_df["run"].astype(str),
        }
    )
    if "calc_mass" in psm_df.columns:
        pin_df["CalcMass"] = psm_df["calc_mass"].astype(float)
    for column in selected_columns:
        pin_df[column] = feature_values[column].to_numpy()
    if not _has_charge_feature(selected_columns):
        pin_df["Charge"] = psm_df["charge"].astype(int)
    pin_df["Peptide"] = [
        to_proforma(sequence, mods, sites)
        for sequence, mods, sites in zip(psm_df["sequence"], psm_df["mods"], psm_df["mod_sites"])
    ]
    pin_df["Proteins"] = psm_df["proteins"].astype(str)
    return pin_df


def mokapot_feature_columns(
    psms: PsmContainer,
    feature_columns: list[str] | tuple[str, ...] | None = None,
) -> tuple[str, ...]:
    """Return the columns Mokapot will use as features."""
    selected_columns = _selected_feature_columns(psms, feature_columns)
    return (
        selected_columns
        if _has_charge_feature(selected_columns)
        else (*selected_columns, "Charge")
    )


def _has_charge_feature(columns: list[str] | tuple[str, ...]) -> bool:
    """Return whether the selected columns already describe precursor charge."""
    return any(
        column == "Charge" or re.fullmatch(r"charge_?\d+(?:_or_more)?", column, re.IGNORECASE)
        for column in columns
    )


def _selected_feature_columns(
    psms: PsmContainer,
    feature_columns: list[str] | tuple[str, ...] | None,
) -> tuple[str, ...]:
    """Return the requested feature columns and always include ``rank``."""
    selected_columns = (
        tuple(feature_columns) if feature_columns is not None else psms.feature_columns
    )
    return selected_columns if "rank" in selected_columns else (*selected_columns, "rank")


def write_pin(
    psms: PsmContainer,
    path: str | Path,
    feature_columns: list[str] | tuple[str, ...] | None = None,
) -> pd.DataFrame:
    """Write PSMs to a PIN file and return the written DataFrame."""
    pin_df = to_mokapot_dataframe(psms, feature_columns=feature_columns)
    pin_df.to_csv(path, sep="\t", index=False)
    return pin_df


def convert_to_mokapot_dataset(
    psms: PsmContainer,
    rescoring_features: list[str] | tuple[str, ...] | None = None,
):
    """Create a Mokapot dataset from the PIN-format PSM DataFrame."""
    selected_columns = (
        tuple(rescoring_features) if rescoring_features is not None else psms.feature_columns
    )
    pin_df = to_mokapot_dataframe(psms, feature_columns=selected_columns)
    dataset = mokapot_lib.read_pin(
        pin_df,
        filename_column="filename",
        calcmass_column="CalcMass" if "CalcMass" in pin_df else None,
        charge_column="Charge" if "Charge" in pin_df else None,
    )
    return dataset


def rescore(
    psms: PsmContainer,
    model=None,
    rescoring_features: list[str] | tuple[str, ...] | None = None,
    test_fdr: float = 0.01,
    rng: int = 1,
    **kwargs,
):
    """Run Mokapot rescoring with the requested features and random seed."""
    dataset = convert_to_mokapot_dataset(psms, rescoring_features=rescoring_features)
    logger.info("Rescoring PSMs with mokapot.")
    return mokapot_lib.brew(
        dataset,
        model=model,
        test_fdr=test_fdr,
        rng=rng,
        **kwargs,
    )
