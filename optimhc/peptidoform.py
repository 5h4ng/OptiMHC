"""Consumer projections of normalized peptidoforms."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import pandas as pd
from pyteomics.proforma import to_proforma as render_proforma


def to_proforma(sequence: object, mods: object, mod_sites: object) -> str:
    """Render normalized modification names as Unimod-accessioned ProForma."""
    sequence = str(sequence)
    residue_tags: dict[int, list[str]] = {}
    nterm: list[str] = []
    cterm: list[str] = []
    if mods:
        names = str(mods).split(";")
        sites = [int(site) for site in str(mod_sites).split(";")]
        if len(names) != len(sites):
            raise ValueError("Modification names and sites must be aligned.")
        identifiers = modification_ids()
        for name, site in zip(names, sites):
            try:
                tag = f"UNIMOD:{identifiers[name]}"
            except KeyError as error:
                raise ValueError(f"Unknown modification '{name}'.") from error
            if site == 0:
                nterm.append(tag)
            elif site == -1:
                cterm.append(tag)
            else:
                residue_tags.setdefault(site, []).append(tag)

    annotated = [
        (residue, residue_tags.get(position, []))
        for position, residue in enumerate(sequence, start=1)
    ]
    return render_proforma(annotated, n_term=nterm or None, c_term=cterm or None)


@lru_cache(maxsize=1)
def modification_ids() -> dict[str, int]:
    table = Path(__file__).parent / "constants" / "modification.tsv"
    modifications = pd.read_csv(table, sep="\t", usecols=["mod_name", "unimod_id"])
    return dict(zip(modifications["mod_name"], modifications["unimod_id"].astype(int)))
