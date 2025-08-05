import logging
from typing import List, Optional, Tuple, Union
import pandas as pd
from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def read_pin(
    pin_files: Union[str, List[str], pd.DataFrame],
    retention_time_column: Optional[str] = None,
    remove_pre_nxt_aa: bool = False,
) -> PsmContainer:
    """
    Read PSMs from a Percolator INput (PIN) file.

    Parameters
    ----------
    pin_files : Union[str, List[str], pd.DataFrame]
        The file path to the PIN file, a list of file paths, or a DataFrame
        containing PIN data.
    retention_time_column : Optional[str], optional
        The column containing the retention time. If None, no retention time
        will be included.

    Returns
    -------
    PsmContainer
        A PsmContainer object containing the PSM data.

    Notes
    -----
    This function:
    1. Reads PIN file(s) into a DataFrame
    2. Identifies required columns (case-insensitive)
    3. Processes scan IDs and hit ranks
    4. Converts data types appropriately
    5. Creates a PsmContainer with the processed data
    """
    logger.info("Reading PIN file(s) into PsmContainer.")
    if isinstance(pin_files, str):
        pin_files = [pin_files]

    pin_df = pd.concat([_read_single_pin_as_df(pin_file) for pin_file in pin_files])
    logger.info(f"Read {len(pin_df)} PSMs from {len(pin_files)} PIN files.")
    logger.debug(pin_df.head())
    logger.debug(pin_df.columns)
    logger.debug(pin_df.iloc[0])

    def find_required_columns(col: str, columns: List[str]) -> str:
        """
        Case-insensitive search for a column in the DataFrame.
        Returns the matching column name with original casing.
        """
        col_lower = col.lower()
        column_map = {c.lower(): c for c in columns}
        if col_lower not in column_map:
            raise ValueError(
                f"Column '{col}' not found in PSM data (case-insensitive)."
            )
        return column_map[col_lower]

    # non-feature columns (case-insensitive search)
    label = find_required_columns("Label", pin_df.columns)
    scan = find_required_columns("ScanNr", pin_df.columns)
    specid = find_required_columns("SpecId", pin_df.columns)
    peptide = find_required_columns("Peptide", pin_df.columns)
    protein = find_required_columns("Proteins", pin_df.columns)

    # Comet: P2PI20160713_pilling_C1RA2_BB72_P1_31_3_1
    # Fragpipe: P2PI20160713_pilling_C1RA2_BB72_P1.3104.3104.2_1

    def parse_specid(specid: str) -> Tuple[str, int]:
        if "_" in specid:
            parts = specid.rsplit("_", 1)
            if len(parts) != 2:
                raise ValueError(f"SpecId format invalid: {specid}")
            unique_id = parts[0]
            hit_rank = int(parts[1])
            return unique_id, hit_rank
        else:
            return specid, 1

    hit_rank = "rank"
    if "rank" in [c.lower() for c in pin_df.columns]:
        pass
    else:
        pin_df[specid], pin_df["rank"] = zip(*pin_df[specid].apply(parse_specid))

    retention_time_column = (
        find_required_columns(retention_time_column, pin_df.columns)
        if retention_time_column
        else None
    )

    # feature columns: columns that are not non-feature columns
    non_feature_columns = [label, scan, specid, peptide, protein]
    feature_columns = [col for col in pin_df.columns if col not in non_feature_columns]

    logger.info(
        f"Columns: label={label}, scan={scan}, specid={specid}, peptide={peptide}, \
                protein={protein}, hit_rank={hit_rank}, retention_time={retention_time_column}, features={feature_columns}"
    )

    pin_df[scan] = pin_df[scan].astype(str)
    pin_df[specid] = pin_df[specid].astype(str)
    pin_df[peptide] = pin_df[peptide].astype(str)
    pin_df[protein] = pin_df[protein].astype(str)
    pin_df[hit_rank] = pin_df[hit_rank].astype(float).astype(int)
    if retention_time_column:
        pin_df[retention_time_column] = pin_df[retention_time_column].astype(float)
    for col in feature_columns:
        pin_df[col] = pin_df[col].astype(float)

    # label = 1 for target, -1 for decoy. Convert to Boolean.
    pin_df[label] = pin_df[label] == "1"
    rescoring_features = {"Original": feature_columns}

    return PsmContainer(
        psms=pin_df,
        label_column=label,
        scan_column=scan,
        spectrum_column=specid,
        ms_data_file_column=None,
        peptide_column=peptide,
        protein_column=protein,
        rescoring_features=rescoring_features,
        hit_rank_column=hit_rank,
        retention_time_column=retention_time_column,
    )


def _read_single_pin_as_df(pin_file: str) -> pd.DataFrame:
    """
    Read a single PIN file into a DataFrame.

    Parameters
    ----------
    pin_file : str
        The file path to the PIN file.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the PSM data.

    Notes
    -----
    This function:
    1. Reads the PIN file header
    2. Processes the proteins column as a tab-separated list
    3. Creates a DataFrame with the processed data
    """
    logger.info(f"Reading PIN file: {pin_file}")
    with open(pin_file, "r") as f:
        header = f.readline().strip().split("\t")
        header_len = len(header)
        data = []
        for line in f:
            parts = line.strip().split("\t")
            proteins_column_num = len(parts) - header_len + 1
            proteins = "\t".join(parts[-proteins_column_num:])
            data.append(parts[: len(parts) - proteins_column_num] + [proteins])
    df = pd.DataFrame(data, columns=header)
    logger.debug(f"Header: {header}")
    return df
