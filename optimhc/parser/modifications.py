"""Match search-engine modifications to the bundled AlphaBase table."""

from __future__ import annotations

import logging
from functools import lru_cache
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


@lru_cache(maxsize=1)
def modification_table() -> pd.DataFrame:
    """Load supported modification names, masses, and Unimod IDs."""
    path = Path(__file__).parents[1] / "constants" / "modification.tsv"
    return pd.read_csv(
        path,
        sep="\t",
        usecols=["mod_name", "unimod_mass", "unimod_id"],
    )


@lru_cache(maxsize=None)
def modification_from_delta(
    delta_mass: float,
    residue: str | None = None,
    site: int | None = None,
    tolerance: float = 0.005,
) -> str:
    """Find a modification name by delta mass and location.

    A missing match raises ``ValueError``. Multiple matches are logged and
    resolved by mass error, terminal scope, then table order.

    Examples
    --------
    >>> modification_from_delta(15.9949, "M")
    'Oxidation@M'
    """
    table = modification_table()
    candidates = table.loc[(table["unimod_mass"] - delta_mass).abs() <= tolerance].copy()
    candidates = _at_location(candidates, residue=residue, site=site)

    if candidates.empty:
        location = "terminal" if site in (0, -1) else residue
        raise ValueError(f"Unknown modification mass {delta_mass:g} at {location}.")
    return _select_best_match(candidates, delta_mass=delta_mass)


@lru_cache(maxsize=None)
def modification_from_unimod(
    unimod_id: int,
    residue: str | None = None,
    site: int | None = None,
) -> str:
    """Find a modification name by Unimod ID and location.

    A missing match raises ``ValueError``. Multiple matches are logged and
    resolved by terminal scope, then table order.

    Examples
    --------
    >>> modification_from_unimod(35, "M")
    'Oxidation@M'
    """
    table = modification_table()
    candidates = table.loc[table["unimod_id"].eq(unimod_id)].copy()
    candidates = _at_location(candidates, residue=residue, site=site)

    if candidates.empty:
        location = "terminal" if site in (0, -1) else residue
        raise ValueError(f"Unknown Unimod accession {unimod_id} at {location}.")
    return _select_best_match(candidates)


def _select_best_match(
    candidates: pd.DataFrame,
    delta_mass: float | None = None,
) -> str:
    """Select a modification using priority rules.

    Priority:
    1. Smallest delta-mass error.
    2. ``Any_[NC]-term`` over ``Protein_[NC]-term``.
    3. Original table order.

    Examples
    --------
    >>> candidates = pd.DataFrame(
    ...     {
    ...         "mod_name": ["Acetyl@Protein_N-term", "Acetyl@Any_N-term"],
    ...         "unimod_mass": [42.010565, 42.010565],
    ...     }
    ... )
    >>> _select_best_match(candidates, 42.010565)
    'Acetyl@Any_N-term'
    """
    candidates = candidates.copy()
    locations = candidates["mod_name"].str.rsplit("@", n=1).str[-1]

    if delta_mass is not None:
        candidates["_mass_error"] = (candidates["unimod_mass"] - delta_mass).abs()
    else:
        candidates["_mass_error"] = 0.0

    any_terminus = locations.str.contains(r"Any_[NC]-term", regex=True, na=False)
    protein_terminus = locations.str.contains(r"Protein_[NC]-term", regex=True, na=False)
    candidates["_terminus_priority"] = 2
    candidates.loc[protein_terminus, "_terminus_priority"] = 1
    candidates.loc[any_terminus, "_terminus_priority"] = 0

    candidates = candidates.sort_values(
        ["_mass_error", "_terminus_priority"],
        kind="stable",
    )
    names = tuple(candidates["mod_name"].astype(str).drop_duplicates())
    selected = names[0]

    if len(names) > 1:
        _warn_ambiguous_mapping(names)
    return selected


@lru_cache(maxsize=None)
def _warn_ambiguous_mapping(names: tuple[str, ...]) -> None:
    """Log each ambiguous set of modification definitions once."""
    logger.warning(
        "Multiple modification definitions match: %s. Using '%s'.",
        ", ".join(names),
        names[0],
    )


def _at_location(
    candidates: pd.DataFrame,
    residue: str | None,
    site: int | None,
) -> pd.DataFrame:
    """Keep table rows that match the requested residue or terminus.

    Examples
    --------
    >>> candidates = pd.DataFrame(
    ...     {
    ...         "mod_name": [
    ...             "Gln->pyro-Glu@Q^Any_N-term",
    ...             "Pyro-carbamidomethyl@C^Any_N-term",
    ...         ]
    ...     }
    ... )
    >>> _at_location(candidates, residue="Q", site=0)["mod_name"].tolist()
    ['Gln->pyro-Glu@Q^Any_N-term']
    """

    locations = candidates["mod_name"].str.rsplit("@", n=1).str[-1]
    if site == 0:
        matches = locations.str.contains("N-term", regex=False)
    elif site == -1:
        matches = locations.str.contains("C-term", regex=False)
    elif residue is not None:
        return candidates[locations.eq(residue)]
    else:
        return candidates

    if residue is not None:
        required_residue = locations.str.extract(r"^([A-Z])\^", expand=False)
        residue_matches = required_residue.isna() | required_residue.eq(residue)
        matches = matches & residue_matches
    return candidates[matches]
