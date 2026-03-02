import logging

import numpy as np
import pandas as pd
from pyteomics import mzml

logger = logging.getLogger(__name__)


def extract_mzml_data(mzml_filename, scan_ids=None):
    """
    Extract scan data from an mzML file.

    Parameters
    ----------
    mzml_filename : str
        The path to the mzML file.
    scan_ids : list[int] or None, optional
        A list of scan IDs to extract. If None, extracts all scans.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the extracted scan data with columns:
        - source: The source file name
        - scan: The scan ID
        - mz: The m/z values array
        - intensity: The intensity values array
        - charge: The charge state
        - retention_time: The retention time

    Notes
    -----
    This function:
    1. Reads the mzML file using pyteomics
    2. Extracts scan data including retention time, charge state, m/z values, and intensities
    3. Filters scans based on provided scan IDs if specified
    4. Returns a DataFrame with the extracted data
    """
    filename = mzml_filename.split("/")[-1].replace(".mzML", "")
    logger.info(f"Extracting scans from {mzml_filename}")

    scan_ids = set(scan_ids) if scan_ids is not None else None

    (
        extracted_scan_ids,
        mzml_filenames,
        intensities,
        mz_values,
        charges,
        retention_times,
    ) = ([], [], [], [], [], [])

    try:
        with mzml.read(mzml_filename) as reader:
            for spectrum in reader:
                try:
                    scan_id = int(spectrum["id"].split("scan=")[-1])

                    if scan_ids is not None and scan_id not in scan_ids:
                        continue

                    mz_array = np.array(spectrum.get("m/z array", []))
                    intensity_array = np.array(spectrum.get("intensity array", []))

                    charge = None
                    try:
                        charge = int(
                            spectrum["precursorList"]["precursor"][0]["selectedIonList"][
                                "selectedIon"
                            ][0]["charge state"]
                        )
                    except (KeyError, ValueError, IndexError):
                        pass

                    retention_time = None
                    try:
                        retention_time = float(spectrum["scanList"]["scan"][0]["scan start time"])
                    except (KeyError, ValueError, IndexError):
                        pass

                    extracted_scan_ids.append(scan_id)
                    mzml_filenames.append(filename)
                    intensities.append(intensity_array)
                    mz_values.append(mz_array)
                    charges.append(charge)
                    retention_times.append(retention_time)

                except Exception as e:
                    logger.warning(f"Skipping scan {scan_id} due to error: {e}")

    except Exception as e:
        logger.error(f"Failed to parse mzML file {mzml_filename}: {e}")
        raise RuntimeError(f"Error processing mzML file {mzml_filename}: {e}")

    data_dict = {
        "source": mzml_filenames,
        "scan": extracted_scan_ids,
        "mz": mz_values,
        "intensity": intensities,
        "charge": charges,
        "retention_time": retention_times,
    }

    scans_df = pd.DataFrame(data_dict)
    scans_df = scans_df.drop_duplicates(subset=["source", "scan"])

    logger.info(f"Successfully extracted {len(scans_df)} scans from {mzml_filename}")

    return scans_df
