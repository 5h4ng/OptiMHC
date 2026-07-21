"""FlashLFQ output adapter for canonical PSMs and Mokapot results."""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

import pandas as pd
from pyteomics.mass import calculate_mass

from optimhc.peptidoform import to_proforma
from optimhc.psm_container import PsmContainer
from optimhc.rescore.mokapot import mokapot_spec_id

logger = logging.getLogger(__name__)

FLASHLFQ_COLUMNS = (
    "File Name",
    "Base Sequence",
    "Full Sequence",
    "Peptide Monoisotopic Mass",
    "Scan Retention Time",
    "Precursor Charge",
    "Protein Accession",
)


def write_flashlfq(
    psms: PsmContainer,
    results,
    path: str | Path,
    *,
    fdr: float,
) -> Path:
    """Write accepted Mokapot peptides in FlashLFQ's tabular format."""
    tables = [
        format_flashlfq(psms, confidence.peptides, fdr=fdr)
        for confidence in _confidence_results(results)
    ]
    output = pd.concat(tables, ignore_index=True) if tables else _empty_output()
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)
    logger.info("FlashLFQ output saved to %s (%d peptides).", output_path, len(output))
    return output_path


def format_flashlfq(
    psms: PsmContainer,
    peptide_results: pd.DataFrame,
    *,
    fdr: float,
) -> pd.DataFrame:
    """Map one Mokapot peptide-confidence table to FlashLFQ columns."""
    if "retention_time" not in psms.df:
        raise ValueError("FlashLFQ export requires canonical column 'retention_time'.")

    required_results = {"SpecId", "Peptide", "mokapot q-value"}
    missing_results = sorted(required_results.difference(peptide_results.columns))
    if missing_results:
        raise ValueError(f"Mokapot peptide results are missing columns: {missing_results}")

    frame = psms.df
    full_sequences = [
        to_proforma(sequence, mods, sites)
        for sequence, mods, sites in zip(frame["sequence"], frame["mods"], frame["mod_sites"])
    ]
    metadata = pd.DataFrame(
        {
            "SpecId": [
                mokapot_spec_id(run, scan, charge)
                for run, scan, charge in zip(frame["run"], frame["scan"], frame["charge"])
            ],
            "Peptide": full_sequences,
            "File Name": frame["run"].map(lambda value: Path(str(value)).name),
            "Base Sequence": frame["sequence"].astype(str),
            "Peptide Monoisotopic Mass": _calculated_masses(frame),
            "Scan Retention Time": pd.to_numeric(
                frame["retention_time"], errors="raise"
            ).astype(float)
            / 60,
            "Precursor Charge": frame["charge"].astype(int),
            "Protein Accession": frame["proteins"].astype(str),
        }
    )
    keys = ["SpecId", "Peptide"]
    metadata = metadata.drop_duplicates()

    accepted = peptide_results.loc[
        pd.to_numeric(peptide_results["mokapot q-value"], errors="raise").le(fdr),
        keys,
    ]
    joined = accepted.merge(
        metadata,
        on=keys,
        how="left",
        sort=False,
        validate="many_to_one",
        indicator=True,
    )
    if not joined["_merge"].eq("both").all():
        missing = joined.loc[joined["_merge"].ne("both"), keys].to_dict("records")
        raise ValueError(f"FlashLFQ could not map accepted peptides to canonical PSMs: {missing}")

    joined["Full Sequence"] = joined["Peptide"]
    return joined.loc[:, FLASHLFQ_COLUMNS]


def _confidence_results(results) -> list:
    if hasattr(results, "peptides"):
        return [results]
    try:
        confidence = list(results)
    except TypeError as error:
        raise ValueError("Mokapot results must expose a peptide confidence table.") from error
    if not all(hasattr(item, "peptides") for item in confidence):
        raise ValueError("Mokapot results must expose a peptide confidence table.")
    return confidence


def _calculated_masses(frame: pd.DataFrame) -> pd.Series:
    calculated = (
        pd.to_numeric(frame["calc_mass"], errors="raise").astype(float)
        if "calc_mass" in frame
        else pd.Series(float("nan"), index=frame.index)
    )
    missing = calculated.isna()
    if missing.any():
        calculated.loc[missing] = [
            _peptide_mass(sequence, mods)
            for sequence, mods in zip(
                frame.loc[missing, "sequence"], frame.loc[missing, "mods"]
            )
        ]
    return calculated


def _peptide_mass(sequence: object, mods: object) -> float:
    mass = float(calculate_mass(sequence=str(sequence)))
    if mods:
        masses = _modification_masses()
        for modification in str(mods).split(";"):
            try:
                mass += masses[modification]
            except KeyError as error:
                raise ValueError(
                    f"Cannot calculate FlashLFQ mass for modification '{modification}'."
                ) from error
    return mass


@lru_cache(maxsize=1)
def _modification_masses() -> dict[str, float]:
    table = Path(__file__).parents[1] / "constants" / "modification.tsv"
    modifications = pd.read_csv(table, sep="\t", usecols=["mod_name", "unimod_mass"])
    return dict(zip(modifications["mod_name"], modifications["unimod_mass"].astype(float)))


def _empty_output() -> pd.DataFrame:
    return pd.DataFrame(columns=FLASHLFQ_COLUMNS)
