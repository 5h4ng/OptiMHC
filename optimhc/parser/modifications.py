"""Normalize search-engine modification masses with AlphaBase definitions."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd


@lru_cache(maxsize=1)
def modification_table() -> pd.DataFrame:
    path = Path(__file__).parents[1] / "constants" / "modification.tsv"
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["mod_name", "unimod_mass", "unimod_id"],
    )


def modification_from_delta(
    delta_mass: float,
    *,
    residue: str | None = None,
    site: int | None = None,
    tolerance: float = 0.01,
) -> str:
    """Resolve one delta mass to an AlphaBase modification name."""
    table = modification_table()
    candidates = table.loc[(table["unimod_mass"] - delta_mass).abs() <= tolerance].copy()
    candidates = _at_location(candidates, residue=residue, site=site)

    if candidates.empty:
        location = "terminal" if site in (0, -1) else residue
        raise ValueError(f"Unknown modification mass {delta_mass:g} at {location}.")
    return str(candidates.iloc[0]["mod_name"])


def modification_from_unimod(
    unimod_id: int,
    *,
    residue: str | None = None,
    site: int | None = None,
) -> str:
    """Resolve a Unimod accession to the normalized AlphaBase name at one site."""
    table = modification_table()
    candidates = table.loc[table["unimod_id"].eq(unimod_id)].copy()
    candidates = _at_location(candidates, residue=residue, site=site)

    if candidates.empty:
        location = "terminal" if site in (0, -1) else residue
        raise ValueError(f"Unknown Unimod accession {unimod_id} at {location}.")
    return str(candidates.iloc[0]["mod_name"])


def _at_location(
    candidates: pd.DataFrame,
    *,
    residue: str | None,
    site: int | None,
) -> pd.DataFrame:
    """Restrict modification definitions to the supported residue or terminus."""

    if site == 0:
        candidates = candidates[candidates["mod_name"].str.contains("N-term", regex=False)]
    elif site == -1:
        candidates = candidates[candidates["mod_name"].str.contains("C-term", regex=False)]
    elif residue is not None:
        candidates = candidates[
            candidates["mod_name"].str.contains(f"@{residue}", regex=False)
            & ~candidates["mod_name"].str.contains("term", case=False, regex=False)
        ]
    return candidates
