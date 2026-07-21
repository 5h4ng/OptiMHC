"""Mokapot boundary for canonical PSM candidates."""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Sequence

import mokapot as mokapot_lib
import numpy as np
import pandas as pd

from optimhc.peptidoform import to_proforma
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def mokapot_spec_id(run: object, scan: object, charge: object) -> str:
    """Render the stable SpecId shared by Mokapot projections and output adapters."""
    return f"{run}.{int(scan)}.{int(scan)}.{int(charge)}"


def to_mokapot_dataframe(
    psms: PsmContainer,
    feature_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Project canonical candidates into one PIN-compatible DataFrame."""
    selected = _selected_feature_columns(psms, feature_columns)
    unknown = [column for column in selected if column not in psms.feature_columns]
    if unknown:
        raise ValueError(f"Unknown rescoring features: {unknown}")

    frame = psms.df
    try:
        numeric_features = frame.loc[:, selected].apply(pd.to_numeric, errors="raise")
    except (TypeError, ValueError) as error:
        raise ValueError("Rescoring features must contain only finite numeric values.") from error
    if not np.isfinite(numeric_features.to_numpy(dtype=float)).all():
        raise ValueError("Rescoring features must contain only finite numeric values.")
    projected = pd.DataFrame(
        {
            "SpecId": [
                mokapot_spec_id(run, scan, charge)
                for run, scan, charge in zip(frame["run"], frame["scan"], frame["charge"])
            ],
            "Label": frame["is_decoy"].map({False: 1, True: -1}).astype(int),
            "ScanNr": frame["scan"].astype(int),
            "filename": frame["run"].astype(str),
        }
    )
    if "calc_mass" in frame.columns:
        projected["CalcMass"] = frame["calc_mass"].astype(float)
    for column in selected:
        projected[column] = numeric_features[column].to_numpy()
    if not _has_charge_feature(selected):
        projected["Charge"] = frame["charge"].astype(int)
    projected["Peptide"] = [
        to_proforma(sequence, mods, sites)
        for sequence, mods, sites in zip(frame["sequence"], frame["mods"], frame["mod_sites"])
    ]
    projected["Proteins"] = frame["proteins"].astype(str)
    return projected


def mokapot_feature_columns(
    psms: PsmContainer,
    feature_columns: Sequence[str] | None = None,
) -> tuple[str, ...]:
    """Return the effective features inferred from the shared PIN projection."""
    selected = _selected_feature_columns(psms, feature_columns)
    return selected if _has_charge_feature(selected) else (*selected, "Charge")


def _has_charge_feature(columns: Sequence[str]) -> bool:
    return any(
        column == "Charge"
        or re.fullmatch(r"charge_?\d+(?:_or_more)?", column, re.IGNORECASE)
        for column in columns
    )


def _selected_feature_columns(
    psms: PsmContainer,
    feature_columns: Sequence[str] | None,
) -> tuple[str, ...]:
    selected = tuple(feature_columns) if feature_columns is not None else psms.feature_columns
    return selected if "rank" in selected else (*selected, "rank")


def write_pin(
    psms: PsmContainer,
    path: str | Path,
    feature_columns: Sequence[str] | None = None,
) -> pd.DataFrame:
    """Write the shared Mokapot projection as a PIN file."""
    projected = to_mokapot_dataframe(psms, feature_columns=feature_columns)
    projected.to_csv(path, sep="\t", index=False)
    return projected


def convert_to_mokapot_dataset(
    psms: PsmContainer,
    rescoring_features: Sequence[str] | None = None,
):
    """Load the shared projection through Mokapot's PIN reader."""
    selected = tuple(rescoring_features) if rescoring_features is not None else psms.feature_columns
    projected = to_mokapot_dataframe(psms, feature_columns=selected)
    dataset = mokapot_lib.read_pin(
        projected,
        filename_column="filename",
        calcmass_column="CalcMass" if "CalcMass" in projected else None,
        charge_column="Charge" if "Charge" in projected else None,
    )
    return dataset


def rescore(
    psms: PsmContainer,
    model=None,
    rescoring_features: Sequence[str] | None = None,
    test_fdr: float = 0.01,
    rng: int = 1,
    **kwargs,
):
    """Rescore canonical candidates deterministically with Mokapot."""
    dataset = convert_to_mokapot_dataset(psms, rescoring_features=rescoring_features)
    logger.info("Rescoring PSMs with mokapot.")
    return mokapot_lib.brew(
        dataset,
        model=model,
        test_fdr=test_fdr,
        rng=rng,
        **kwargs,
    )
