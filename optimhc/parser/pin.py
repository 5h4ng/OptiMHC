"""Reader for Percolator input (PIN) files."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from optimhc.parser.modifications import modification_from_delta, modification_from_unimod
from optimhc.psm_container import PsmContainer
from optimhc.utils import strip_flanking_and_charge


def read_pin(
    pin_file: str | Path,
    retention_time_column: str | None = None,
) -> PsmContainer:
    """Read one PIN file into the canonical PSM table."""
    path = Path(pin_file)
    pin = _read_single_pin_as_df(path)
    columns = {column.lower(): column for column in pin.columns}

    spec_id = _required_column(columns, "specid")
    label = _required_column(columns, "label")
    scan = _required_column(columns, "scannr")
    peptide = _required_column(columns, "peptide")
    proteins = _required_column(columns, "proteins")
    rank = columns.get("rank")
    filename = columns.get("filename")
    exp_mass = columns.get("expmass")
    calc_mass = columns.get("calcmass")

    if retention_time_column is not None:
        rt = _required_column(columns, retention_time_column.lower())
    else:
        rt = next(
            (
                columns[name]
                for name in ("retention_time", "retentiontime", "ret_time")
                if name in columns
            ),
            None,
        )

    scan_values = pd.to_numeric(pin[scan], errors="raise").astype(int)
    rank_values = (
        pd.to_numeric(pin[rank], errors="raise").astype(int)
        if rank is not None
        else pin[spec_id].map(_rank_from_spec_id).astype(int)
    )
    charge_values = _read_charge(pin)
    fallback_run = _normalize_run(path.stem)
    run_values = [
        _normalize_run(value) if filename is not None else _run_from_spec_id(sid, sn, fallback_run)
        for sid, sn, value in zip(
            pin[spec_id],
            scan_values,
            pin[filename] if filename is not None else [fallback_run] * len(pin),
        )
    ]
    peptidoforms = pin[peptide].map(_parse_pin_peptide)
    labels = pd.to_numeric(pin[label], errors="raise")
    if not labels.isin((-1, 1)).all():
        raise ValueError("PIN labels must be exactly -1 or 1.")

    metadata = {
        spec_id,
        label,
        scan,
        peptide,
        proteins,
        *(column for column in (rank, filename, exp_mass, calc_mass, rt) if column is not None),
    }
    feature_columns = tuple(column for column in pin.columns if column not in metadata)
    if rank is None or rank not in feature_columns:
        feature_columns = (*feature_columns, "rank")

    canonical = pd.DataFrame(
        {
            "psm_id": range(len(pin)),
            "run": run_values,
            "scan": scan_values.to_numpy(),
            "rank": rank_values.to_numpy(),
            "sequence": [item[0] for item in peptidoforms],
            "mods": [item[1] for item in peptidoforms],
            "mod_sites": [item[2] for item in peptidoforms],
            "charge": charge_values.to_numpy(),
            "proteins": pin[proteins].map(_normalize_proteins),
            "is_decoy": labels.eq(-1).to_numpy(),
        }
    )
    if rt is not None:
        retention_time = pd.to_numeric(pin[rt], errors="raise")
        if rt.lower() == "retentiontime":
            retention_time = retention_time * 60
        canonical["retention_time"] = retention_time
    if exp_mass is not None:
        canonical["exp_mass"] = pd.to_numeric(pin[exp_mass], errors="raise")
    if calc_mass is not None:
        canonical["calc_mass"] = pd.to_numeric(pin[calc_mass], errors="raise")
    for column in feature_columns:
        if column == "rank":
            continue
        canonical[column] = pd.to_numeric(pin[column], errors="raise").to_numpy()

    return PsmContainer(canonical, feature_columns=feature_columns)


def _read_single_pin_as_df(pin_file: str | Path) -> pd.DataFrame:
    with Path(pin_file).open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line in handle:
            values = line.rstrip("\n").split("\t")
            protein_count = len(values) - len(header) + 1
            rows.append(values[: len(values) - protein_count] + ["\t".join(values[-protein_count:])])
    return pd.DataFrame(rows, columns=header)


def _required_column(columns: dict[str, str], name: str) -> str:
    try:
        return columns[name]
    except KeyError as error:
        raise ValueError(f"Column '{name}' not found in PIN data.") from error


def _read_charge(pin: pd.DataFrame) -> pd.Series:
    direct = next((column for column in pin.columns if column.lower() == "charge"), None)
    if direct is not None:
        return pd.to_numeric(pin[direct], errors="raise").astype(int)

    charge_columns = {
        column: int(match.group(1))
        for column in pin.columns
        if (match := re.fullmatch(r"charge_?(\d+)(?:_or_more)?", column, re.IGNORECASE))
    }
    if not charge_columns:
        raise ValueError("PIN data must contain charge or one-hot charge columns.")

    def selected_charge(row: pd.Series) -> int:
        selected = [charge for column, charge in charge_columns.items() if float(row[column]) == 1]
        if len(selected) != 1:
            raise ValueError("Each PIN row must select exactly one charge state.")
        return selected[0]

    return pin.apply(selected_charge, axis=1).astype(int)


def _rank_from_spec_id(spec_id: object) -> int:
    match = re.search(r"_(\d+)$", str(spec_id))
    return int(match.group(1)) if match else 1


def _run_from_spec_id(spec_id: object, scan: int, fallback: str) -> str:
    text = str(spec_id)
    marker = f".{scan}."
    return _normalize_run(text.split(marker, 1)[0]) if marker in text else fallback


def _normalize_run(value: object) -> str:
    name = Path(str(value)).name
    for suffix in (".mzML", ".mzml", ".raw", ".RAW"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _parse_pin_peptide(value: object) -> tuple[str, str, str]:
    peptide = strip_flanking_and_charge(str(value))
    sequence: list[str] = []
    mods: list[str] = []
    sites: list[str] = []

    nterm = re.match(r"^\[UNIMOD:(\d+)\]-", peptide, re.IGNORECASE)
    if nterm:
        mods.append(modification_from_unimod(int(nterm.group(1)), site=0))
        sites.append("0")
        peptide = peptide[nterm.end() :]
    cterm = re.search(r"-\[UNIMOD:(\d+)\]$", peptide, re.IGNORECASE)
    cterm_mod = None
    if cterm:
        cterm_mod = modification_from_unimod(int(cterm.group(1)), site=-1)
        peptide = peptide[: cterm.start()]

    index = 0
    while index < len(peptide):
        terminal = None
        if peptide.startswith("n[", index) and index == 0:
            terminal = 0
        elif peptide.startswith("c[", index) and sequence:
            terminal = -1
        if terminal is not None:
            end = peptide.find("]", index + 2)
            if end < 0:
                raise ValueError(f"Invalid modified peptide: {value}")
            try:
                mass = float(peptide[index + 2 : end])
            except ValueError as error:
                raise ValueError(f"Unknown terminal modification in '{value}'.") from error
            mods.append(modification_from_delta(mass, site=terminal))
            sites.append(str(terminal))
            index = end + 1
            continue
        if peptide[index] == "[":
            end = peptide.find("]", index)
            if end < 0 or not sequence:
                raise ValueError(f"Invalid modified peptide: {value}")
            annotation = peptide[index + 1 : end]
            unimod = re.fullmatch(r"UNIMOD:(\d+)", annotation, re.IGNORECASE)
            if unimod:
                mod_name = modification_from_unimod(
                    int(unimod.group(1)), residue=sequence[-1]
                )
                mod_site = len(sequence)
            else:
                try:
                    mass = float(annotation)
                except ValueError as error:
                    raise ValueError(f"Unknown modification '{annotation}'.") from error
                try:
                    mod_name = modification_from_delta(mass, residue=sequence[-1])
                    mod_site = len(sequence)
                except ValueError:
                    if len(sequence) != 1:
                        raise
                    mod_name = modification_from_delta(mass, site=0)
                    mod_site = 0
            mods.append(mod_name)
            sites.append(str(mod_site))
            index = end + 1
        else:
            sequence.append(peptide[index])
            index += 1
    if cterm_mod is not None:
        mods.append(cterm_mod)
        sites.append("-1")
    return "".join(sequence), ";".join(mods), ";".join(sites)


def _normalize_proteins(value: object) -> str:
    return ";".join(part for part in re.split(r"[;\t]+", str(value)) if part)
