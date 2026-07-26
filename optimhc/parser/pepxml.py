"""Read PSMs from pepXML files."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from lxml import etree
from pyteomics.mass import std_aa_mass

from optimhc.parser.modifications import modification_from_delta
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)
PROTON = 1.00727646677


def read_pepxml(
    pepxml_file: str | Path,
    decoy_prefix: str = "DECOY_",
) -> PsmContainer:
    """Read every search hit from one pepXML file.

    Each ``spectrum_query`` must provide scan, charge, retention time in seconds,
    and precursor mass. Each ``search_hit`` must provide peptide, protein, and
    calculated mass. Modification masses must match the bundled modification
    table. Percolator and PeptideProphet output is not accepted as search input.
    """
    records = _parse_pepxml(Path(pepxml_file), decoy_prefix)
    psms = pd.DataFrame.from_records(records)
    illegal = {"Percolator q-Value", "Percolator PEP", "Percolator SVMScore"}
    if illegal.intersection(psms.columns):
        raise ValueError("Percolator or PeptideProphet output cannot be rescored as raw pepXML.")

    psms["mass_diff"] = psms["exp_mass"] - psms["calc_mass"]
    exp_mz = psms["exp_mass"] / psms["charge"] + PROTON
    calc_mz = psms["calc_mass"] / psms["charge"] + PROTON
    psms["abs_mz_diff"] = (exp_mz - calc_mz).abs()
    # assert tot_num_ions is never zero for pepxml
    if {"num_matched_ions", "tot_num_ions"}.issubset(psms.columns):
        psms["matched_ions_ratio"] = psms["num_matched_ions"] / psms["tot_num_ions"]
    if "num_matched_peptides" in psms.columns:
        psms["num_matched_peptides"] = np.log10(psms["num_matched_peptides"])

    charge_features = pd.get_dummies(psms["charge"], prefix="charge", dtype=float)
    psms = pd.concat([psms, charge_features], axis=1)
    psm_columns = {
        "run",
        "scan",
        "rank",
        "sequence",
        "mods",
        "mod_sites",
        "charge",
        "proteins",
        "is_decoy",
        "retention_time",
        "exp_mass",
        "calc_mass",
    }
    feature_columns = tuple(column for column in psms.columns if column not in psm_columns)
    for column in feature_columns:
        psms[column] = _log_feature(psms[column])

    psms.insert(0, "psm_id", range(len(psms)))
    ordered = [
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
        "retention_time",
        "exp_mass",
        "calc_mass",
        *feature_columns,
    ]
    return PsmContainer(psms.loc[:, ordered], feature_columns=feature_columns)


def _parse_pepxml(path: Path, decoy_prefix: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    try:
        spectra = etree.iterparse(str(path), events=("end",), tag="{*}spectrum_query")
        for _, spectrum in spectra:
            run_summary = spectrum.getparent()
            run = _normalize_run(run_summary.get("base_name"))
            records.extend(_parse_spectrum(spectrum, run, decoy_prefix))
            # TODO: clear processed elements for large pepXML files.
    except etree.XMLSyntaxError as error:
        raise ValueError(f"{path} is not a valid pepXML file.") from error
    if not records:
        raise ValueError(f"No peptide-spectrum matches found in {path}.")
    return records


def _parse_spectrum(
    spectrum: etree._Element,
    run: str,
    decoy_prefix: str,
) -> list[dict[str, object]]:
    """Return one row for every search hit in a spectrum query."""
    spectrum_fields = {
        "run": run,
        "scan": int(spectrum.get("end_scan")),
        "charge": int(spectrum.get("assumed_charge")),
        "retention_time": float(spectrum.get("retention_time_sec")),
        "exp_mass": float(spectrum.get("precursor_neutral_mass")),
    }
    records = []
    for result in spectrum.iter("{*}search_result"):
        for hit in result.iter("{*}search_hit"):
            records.append(_parse_hit(hit, spectrum_fields, decoy_prefix))
    return records


def _parse_hit(
    hit: etree._Element,
    spectrum_fields: dict[str, object],
    decoy_prefix: str,
) -> dict[str, object]:
    """Convert one ``search_hit`` element to a PSM row.

    A row is a decoy only when every listed protein starts with
    ``decoy_prefix``. Missing ``hit_rank`` values are stored as rank 1.
    """
    sequence = str(hit.get("peptide"))
    proteins = [str(hit.get("protein")).split(" ", 1)[0]]
    proteins.extend(
        str(element.get("protein")).split(" ", 1)[0]
        for element in hit.iter("{*}alternative_protein")
    )
    record = {
        **spectrum_fields,
        "rank": int(hit.get("hit_rank", 1)),
        "sequence": sequence,
        "mods": "",
        "mod_sites": "",
        "calc_mass": float(hit.get("calc_neutral_pep_mass")),
        "proteins": ";".join(proteins),
        "is_decoy": all(protein.startswith(decoy_prefix) for protein in proteins),
    }
    integer_attributes = {
        "num_missed_cleavages": "missed_cleavages",
        "num_tol_term": "ntt",
        "num_matched_peptides": "num_matched_peptides",
        "num_matched_ions": "num_matched_ions",
        "tot_num_ions": "tot_num_ions",
    }
    for attribute, column in integer_attributes.items():
        if hit.get(attribute) is not None:
            record[column] = int(hit.get(attribute))
    mods, sites = _parse_modifications(hit, sequence)
    record["mods"] = ";".join(mods)
    record["mod_sites"] = ";".join(str(site) for site in sites)
    for score in hit.iter("{*}search_score"):
        record[str(score.get("name"))] = score.get("value")
    return record


def _parse_modifications(hit: etree._Element, sequence: str) -> tuple[list[str], list[int]]:
    """Return modification names and one-based residue or terminal sites.

    Sites use ``0`` for the N-terminus and ``-1`` for the C-terminus. Reported
    residue masses are converted to delta masses before table lookup.

    Examples
    --------
    >>> hit = etree.fromstring(
    ...     '<search_hit><modification_info>'
    ...     '<mod_aminoacid_mass position="1" mass="147.0354"/>'
    ...     '</modification_info></search_hit>'
    ... )
    >>> _parse_modifications(hit, "M")
    (['Oxidation@M'], [1])
    """
    info = next(hit.iter("{*}modification_info"), None)
    if info is None:
        return [], []

    mods: list[str] = []
    sites: list[int] = []
    if info.get("mod_nterm_mass") is not None:
        delta = float(info.get("mod_nterm_mass")) - 1.00782503223
        mods.append(modification_from_delta(delta, residue=sequence[0], site=0))
        sites.append(0)
    for modification in info.iter("{*}mod_aminoacid_mass"):
        position = int(modification.get("position"))
        residue = sequence[position - 1]
        delta = float(modification.get("mass")) - std_aa_mass[residue]
        # some engines may report n-terminal mods as position = 1 residue mods
        # fall back to site = 0 when residue mod lookup fails at position 1
        try:
            mod_name = modification_from_delta(delta, residue=residue)
            mod_site = position
        except ValueError:
            if position != 1:
                raise
            mod_name = modification_from_delta(delta, residue=residue, site=0)
            mod_site = 0
        mods.append(mod_name)
        sites.append(mod_site)
    if info.get("mod_cterm_mass") is not None:
        delta = float(info.get("mod_cterm_mass")) - 17.00273965163
        mods.append(modification_from_delta(delta, residue=sequence[-1], site=-1))
        sites.append(-1)
    return mods, sites


def _normalize_run(value: object) -> str:
    """Return a run basename without a raw-data file extension."""
    name = Path(str(value)).name
    for suffix in (".mzML", ".mzml"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _log_feature(column: pd.Series) -> pd.Series:
    """Apply the existing p-value/E-value log rule to a numeric feature.

    Values must be numeric. Non-negative columns spanning at least four orders
    of magnitude are log-transformed; other columns are returned as floats.

    Examples
    --------
    >>> _log_feature(pd.Series([1e-6, 1e-1])).tolist()
    [-6.0, -1.0]
    """
    if pd.api.types.is_bool_dtype(column):
        return column.astype(float)
    text = column.astype(str).str.lower()
    numeric = pd.to_numeric(text, errors="raise")
    if text.str.contains("e").any() and numeric.gt(0).all():
        split = text.str.split("e", expand=True)
        root = pd.to_numeric(split.loc[:, 0])
        power = pd.to_numeric(split.loc[:, 1].fillna("0")).astype(int)
        zero = root.eq(0)
        root.loc[zero] = 1
        if zero.any() and (~zero).any():
            power.loc[zero] = power.loc[~zero].min()
        if abs(power.max() - power.min()) >= 4:
            logger.info("Log-transformed feature '%s'.", column.name)
            return np.log10(root) + power
    if numeric.min() >= 0 and numeric.max() > 0:
        positive = numeric[numeric > 0]
        if not positive.empty and numeric.max() / positive.min() >= 10_000:
            logged = np.log10(numeric.clip(lower=positive.min()))
            zero_mask = numeric.eq(0)
            if zero_mask.any():
                logged[zero_mask] = logged[~zero_mask].min() - 1
            return logged
    return numeric.astype(float)
