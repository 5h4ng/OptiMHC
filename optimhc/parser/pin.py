"""Read PSMs from Percolator input (PIN) files."""

from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd

from optimhc.parser.modifications import modification_from_delta, modification_from_unimod
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)
_NUMERIC_MASS = r"[+-]?\d+(?:\.\d+)?"


def read_pin(
    pin_file: str | Path,
    retention_time_column: str | None = None,
) -> PsmContainer:
    """Read one PIN file.

    The file must contain ``SpecId``, ``Label``, ``ScanNr``, ``Peptide``, and
    ``Proteins`` plus either ``Charge`` or supported one-hot charge columns.
    Labels must be ``1`` or ``-1``. Unknown non-numeric columns are omitted
    from the PSM table and reported as warnings. If the largest retention time
    is below 500, the values are treated as minutes and converted to seconds.
    """
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
        rt = None
        for name in ("retention_time", "retentiontime", "ret_time", "rt"):
            if name in columns:
                rt = columns[name]
                break

    scan_values = pd.to_numeric(pin[scan], errors="raise").astype(int)
    if rank is not None:
        rank_values = pd.to_numeric(pin[rank], errors="raise").astype(int)
    else:
        rank_values = pin[spec_id].map(_rank_from_spec_id).astype(int)

    charge_values = _read_charge(pin)
    pin_file_run = _normalize_run(path.stem)

    run_values = []
    for row_number, (spec_id_value, scan_number) in enumerate(zip(pin[spec_id], scan_values)):
        if filename is not None:
            run = _normalize_run(pin[filename].iloc[row_number])
        else:
            # FragPipe SpecId starts with <run>.<scan>.
            spec_id_text = str(spec_id_value)
            scan_marker = f".{scan_number}."
            if scan_marker in spec_id_text:
                run = _normalize_run(spec_id_text.split(scan_marker, 1)[0])
            else:
                run = pin_file_run
        run_values.append(run)

    parsed_peptides = pin[peptide].map(_parse_pin_peptide)
    labels = pd.to_numeric(pin[label], errors="raise")
    if not labels.isin((-1, 1)).all():
        raise ValueError("PIN labels must be exactly -1 or 1.")

    metadata = {spec_id, label, scan, peptide, proteins}
    for column in (rank, filename, exp_mass, calc_mass, rt):
        if column is not None:
            metadata.add(column)
    candidate_feature_columns = tuple(column for column in pin.columns if column not in metadata)

    psm_df = pd.DataFrame(
        {
            "psm_id": range(len(pin)),
            "run": run_values,
            "scan": scan_values.to_numpy(),
            "rank": rank_values.to_numpy(),
            "sequence": [item[0] for item in parsed_peptides],
            "mods": [item[1] for item in parsed_peptides],
            "mod_sites": [item[2] for item in parsed_peptides],
            "charge": charge_values.to_numpy(),
            "proteins": pin[proteins].map(_normalize_proteins),
            "is_decoy": labels.eq(-1).to_numpy(),
        }
    )
    if rt is not None:
        retention_time = pd.to_numeric(pin[rt], errors="raise")
        if retention_time.max() < 500:
            logger.info("Converting PIN retention times from minutes to seconds.")
            retention_time = retention_time * 60
        psm_df["retention_time"] = retention_time
    if exp_mass is not None:
        psm_df["exp_mass"] = pd.to_numeric(pin[exp_mass], errors="raise")
    if calc_mass is not None:
        psm_df["calc_mass"] = pd.to_numeric(pin[calc_mass], errors="raise")
    feature_columns = []
    for column in candidate_feature_columns:
        try:
            values = pd.to_numeric(pin[column], errors="raise")
        except (TypeError, ValueError):
            logger.warning(
                "Ignoring non-numeric PIN column '%s'; it will not be used for rescoring.",
                column,
            )
            continue
        psm_df[column] = values.to_numpy()
        feature_columns.append(column)

    feature_columns.append("rank")

    return PsmContainer(psm_df, feature_columns=tuple(feature_columns))


def _read_single_pin_as_df(pin_file: str | Path) -> pd.DataFrame:
    """Read a tab-separated PIN file whose final field contains proteins.

    Extra tab-separated values after the fixed columns are treated as additional
    protein accessions.
    """
    with Path(pin_file).open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = []
        for line in handle:
            values = line.rstrip("\n").split("\t")
            protein_count = len(values) - len(header) + 1
            rows.append(
                values[: len(values) - protein_count] + ["\t".join(values[-protein_count:])]
            )
    return pd.DataFrame(rows, columns=header)


def _required_column(columns: dict[str, str], name: str) -> str:
    try:
        return columns[name]
    except KeyError as error:
        raise ValueError(f"Column '{name}' not found in PIN data.") from error


def _read_charge(pin: pd.DataFrame) -> pd.Series:
    """Read charge values from one integer column or one-hot columns.

    When one-hot columns are used, each row must select exactly one charge.

    Examples
    --------
    >>> _read_charge(pd.DataFrame({"Charge": [3]})).tolist()
    [3]
    >>> _read_charge(
    ...     pd.DataFrame({"charge_2": [0], "charge_3": [1]})
    ... ).tolist()
    [3]
    """
    charge_column = None
    for column in pin.columns:
        if column.lower() == "charge":
            charge_column = column
            break
    if charge_column is not None:
        return pd.to_numeric(pin[charge_column], errors="raise").astype(int)

    charge_columns = {}
    for column in pin.columns:
        match = re.fullmatch(r"charge_?(\d+)(?:_or_more)?", column, re.IGNORECASE)
        if match:
            charge_columns[column] = int(match.group(1))
    if not charge_columns:
        raise ValueError("PIN data must contain charge or one-hot charge columns.")

    def selected_charge(row: pd.Series) -> int:
        selected_charges = [
            charge for column, charge in charge_columns.items() if float(row[column]) == 1
        ]
        if len(selected_charges) != 1:
            raise ValueError("Each PIN row must select exactly one charge state.")
        return selected_charges[0]

    return pin.apply(selected_charge, axis=1).astype(int)


def _rank_from_spec_id(spec_id: object) -> int:
    """Read a trailing ``_<rank>`` value, or return rank 1 when absent.

    Examples
    --------
    >>> _rank_from_spec_id("run.10.10.2_3")
    3
    >>> _rank_from_spec_id("run.10.10.2")
    1
    """
    match = re.search(r"_(\d+)$", str(spec_id))
    return int(match.group(1)) if match else 1


def _normalize_run(value: object) -> str:
    """Return a run basename without a raw-data file extension.

    Examples
    --------
    >>> _normalize_run("/data/run_a.mzML")
    'run_a'
    """
    name = Path(str(value)).name
    for suffix in (".mzML", ".mzml"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def _strip_pin_flanks_and_charge(value: object) -> str:
    """Remove optional flanks and a trailing charge number.

    Examples
    --------
    >>> _strip_pin_flanks_and_charge("R.Q[-17.0265]PEPTIDE4.S")
    'Q[-17.0265]PEPTIDE'
    """

    peptide = str(value)
    flanked = re.fullmatch(r"[A-Z-]\.(.+)\.[A-Z-]", peptide, re.IGNORECASE)
    if flanked:
        peptide = flanked.group(1)
    if peptide and peptide[-1].isdigit():
        peptide = peptide[:-1]
    return peptide


def _parse_pin_peptide(value: object) -> tuple[str, str, str]:
    """Parse a PIN peptide into sequence, modification names, and sites.

    Flanking residues are optional. External PIN input supports numeric delta
    masses on residues and both common terminal styles:

    - Comet/MSFragger: ``n[mass]PEPTIDE`` and ``PEPTIDEc[mass]``.
    - Sage: ``[mass]-PEPTIDE`` and ``PEPTIDE-[mass]``.

    Numeric masses and ``UNIMOD:<id>`` annotations must match the bundled
    modification table. Symbol modifications are rejected.

    Returns the unmodified sequence and aligned semicolon-separated ``mods``
    and ``mod_sites`` strings. Site ``0`` means the N-terminus, ``-1`` means
    the C-terminus, and ``1..N`` are one-based residue positions. Malformed or
    unknown annotations raise ``ValueError``.

    Examples
    --------
    >>> _parse_pin_peptide("K.n[42.0106]PEPM[15.9949]K.R")
    ('PEPMK', 'Acetyl@Any_N-term;Oxidation@M', '0;4')
    """
    peptide = _strip_pin_flanks_and_charge(value)
    sequence: list[str] = []
    mods: list[str] = []
    sites: list[str] = []

    # Read one N-terminal modification before scanning residues.
    nterm_unimod = re.match(r"^\[UNIMOD:(\d+)\]-", peptide, re.IGNORECASE)
    nterm_mass = re.match(rf"^n\[({_NUMERIC_MASS})\]", peptide)
    if nterm_mass is None:
        nterm_mass = re.match(rf"^\[({_NUMERIC_MASS})\]-", peptide)
    if nterm_unimod:
        peptide = peptide[nterm_unimod.end() :]
        residue = peptide[0] if peptide else None
        mods.append(
            modification_from_unimod(
                int(nterm_unimod.group(1)),
                residue=residue,
                site=0,
            )
        )
        sites.append("0")
    elif nterm_mass:
        peptide = peptide[nterm_mass.end() :]
        residue = peptide[0] if peptide else None
        mods.append(
            modification_from_delta(
                float(nterm_mass.group(1)),
                residue=residue,
                site=0,
            )
        )
        sites.append("0")

    # Remove one C-terminal modification and append it after residue mods.
    cterm_unimod = re.search(r"-\[UNIMOD:(\d+)\]$", peptide, re.IGNORECASE)
    cterm_mass = re.search(rf"c\[({_NUMERIC_MASS})\]$", peptide)
    if cterm_mass is None:
        cterm_mass = re.search(rf"-\[({_NUMERIC_MASS})\]$", peptide)
    cterm_unimod_id = None
    cterm_delta = None
    if cterm_unimod:
        cterm_unimod_id = int(cterm_unimod.group(1))
        peptide = peptide[: cterm_unimod.start()]
    elif cterm_mass:
        cterm_delta = float(cterm_mass.group(1))
        peptide = peptide[: cterm_mass.start()]

    index = 0
    while index < len(peptide):
        # A bracket after a residue contains a Unimod ID or delta mass.
        if peptide[index] == "[":
            end = peptide.find("]", index)
            if end < 0 or not sequence:
                raise ValueError(f"Invalid modified peptide: {value}")
            annotation = peptide[index + 1 : end]
            unimod = re.fullmatch(r"UNIMOD:(\d+)", annotation, re.IGNORECASE)
            if unimod:
                mod_name = modification_from_unimod(int(unimod.group(1)), residue=sequence[-1])
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
                    # Some PIN files attach an N-terminal, residue-specific
                    # mass to residue 1, for example Q[-17.0265].
                    if len(sequence) != 1:
                        raise
                    mod_name = modification_from_delta(
                        mass,
                        residue=sequence[-1],
                        site=0,
                    )
                    mod_site = 0
            mods.append(mod_name)
            sites.append(str(mod_site))
            index = end + 1
        else:
            if re.fullmatch(r"[A-Z]", peptide[index], re.IGNORECASE) is None:
                raise ValueError(f"Unsupported PIN peptide syntax: {value}")
            sequence.append(peptide[index])
            index += 1
    if cterm_unimod_id is not None:
        mods.append(
            modification_from_unimod(
                cterm_unimod_id,
                residue=sequence[-1],
                site=-1,
            )
        )
        sites.append("-1")
    elif cterm_delta is not None:
        mods.append(
            modification_from_delta(
                cterm_delta,
                residue=sequence[-1],
                site=-1,
            )
        )
        sites.append("-1")
    return "".join(sequence), ";".join(mods), ";".join(sites)


def _normalize_proteins(value: object) -> str:
    """Join tab- or semicolon-separated protein names with semicolons.

    Examples
    --------
    >>> _normalize_proteins("P1\tP2")
    'P1;P2'
    """
    return ";".join(part for part in re.split(r"[;\t]+", str(value)) if part)
