"""Convert OptiMHC peptide data to ProForma."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
from pyteomics.proforma import to_proforma as format_proforma


def to_proforma(sequence: object, mods: object, mod_sites: object) -> str:
    """Convert a sequence and its modification columns to ProForma.

    ``mods`` and ``mod_sites`` must be aligned semicolon-separated values.
    Sites use ``0`` for the N-terminus, ``-1`` for the C-terminus, and
    one-based positions for residues. Every modification name must exist in
    the bundled modification table.
    """
    sequence = str(sequence)
    residue_tags: dict[int, list[str]] = {}
    nterm: list[str] = []
    cterm: list[str] = []
    if mods:
        mod_names = str(mods).split(";")
        site_values = [int(site) for site in str(mod_sites).split(";")]
        if len(mod_names) != len(site_values):
            raise ValueError("Modification names and sites must be aligned.")
        unimod_ids = modification_ids()
        for name, site in zip(mod_names, site_values):
            try:
                unimod_tag = f"UNIMOD:{unimod_ids[name]}"
            except KeyError as error:
                raise ValueError(f"Unknown modification '{name}'.") from error
            if site == 0:
                nterm.append(unimod_tag)
            elif site == -1:
                cterm.append(unimod_tag)
            else:
                residue_tags.setdefault(site, []).append(unimod_tag)

    residues_with_mods = [
        (residue, residue_tags.get(position, []))
        for position, residue in enumerate(sequence, start=1)
    ]
    return format_proforma(residues_with_mods, n_term=nterm or None, c_term=cterm or None)


@lru_cache(maxsize=1)
def modification_ids() -> dict[str, int]:
    """Return the Unimod ID for each supported modification name."""
    table = Path(__file__).parent / "constants" / "modification.tsv"
    modifications = pd.read_csv(table, sep="\t", usecols=["mod_name", "unimod_id"])
    return dict(zip(modifications["mod_name"], modifications["unimod_id"].astype(int)))
