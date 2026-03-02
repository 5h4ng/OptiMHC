# rescore/mokapot.py
# A wrapper around mokapot for rescoring PSMs and converting to flashLFQ format

import logging
import os
from pathlib import Path
from typing import List

import mokapot
import pandas as pd
from mokapot import LinearPsmDataset

from optimhc.psm_container import PsmContainer

logger = logging.getLogger(__name__)


def rescore(
    psms: PsmContainer,
    model=None,
    rescoring_features: List[str] = None,
    test_fdr: float = 0.01,
    **kwargs,
):
    """
    Rescore PSMs using mokapot.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing PSM data.
    model : object, optional
        A trained model for rescoring PSMs.
    rescoring_features : List[str], optional
        A list of feature names to use for rescoring.
    test_fdr : float, optional
        The FDR threshold for testing the model. Default is 0.01.
    **kwargs : dict
        Additional keyword arguments for mokapot.brew.

    Returns
    -------
    tuple
        A tuple containing:
        - Confidence object or list of Confidence objects:
          An object or a list of objects containing the confidence estimates at various levels
          (i.e. PSMs, peptides) when assessed using the learned score. If a list, they will be
          in the same order as provided in the psms parameter.
        - list of Model objects:
          The learned Model objects, one for each fold.

    Notes
    -----
    This function:
    1. Converts the PsmContainer to a mokapot dataset
    2. Runs mokapot.brew with the specified parameters
    3. Returns the results and models
    """
    psms = convert_to_mokapot_dataset(psms, rescoring_features=rescoring_features)
    logger.info("Rescoring PSMs with mokapot.")
    results, models = mokapot.brew(psms, model=model, test_fdr=test_fdr, **kwargs)
    return results, models


def convert_to_mokapot_dataset(
    psms: PsmContainer, rescoring_features: List[str] = None
) -> LinearPsmDataset:
    """
    Convert a PsmContainer to a LinearPsmDataset for use with mokapot.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing PSM data.
    rescoring_features : List[str], optional
        A list of feature names to use for rescoring.
        If not provided, uses all features from the PsmContainer.

    Returns
    -------
    LinearPsmDataset
        A LinearPsmDataset object for use with mokapot.

    Raises
    ------
    ValueError
        If any of the specified rescoring features are not found in the PSM data.

    Notes
    -----
    This function:
    1. Extracts all features from the PsmContainer
    2. Validates the specified rescoring features
    3. Creates a LinearPsmDataset with the appropriate columns and data
    """

    # rescoring_features: Dict[str, List[str]] -> Tuple[str]
    feature_columns = [col for features in psms.rescoring_features.values() for col in features]

    if rescoring_features is None:
        rescoring_features = feature_columns
    else:
        for feature in rescoring_features:
            if feature not in feature_columns:
                raise ValueError(f"Feature '{feature}' not found in the PSM data.")

    dataset = LinearPsmDataset(
        psms.psms,
        target_column=psms.label_column,
        spectrum_columns=psms.spectrum_column,
        peptide_column=psms.peptide_column,
        protein_column=psms.protein_column,
        feature_columns=rescoring_features,
        filename_column=psms.ms_data_file_column,
        calcmass_column=psms.calculated_mass_column,
        charge_column=psms.charge_column,
        scan_column=psms.scan_column,
        rt_column=psms.retention_time_column,
    )

    return dataset


# Adapted from mokapot source code
# https://github.com/wfondrie/mokapot


def to_flashLFQ(results, output_dir, file_name):
    logger.info("Saving results in FlashLFQ format.")
    try:
        assert not isinstance(results, str)
        iter(results)
    except TypeError:
        results = [results]
    except AssertionError:
        raise ValueError("'results' should be a Confidence object, not a string.")
    flashlfq = pd.concat([_format_flashlfq(c) for c in results])
    flashlfq.to_csv(os.path.join(output_dir, file_name), sep="\t", index=False)
    return os.path.join(output_dir, file_name)


def _format_flashlfq(conf):
    # Do some error checking for the required columns:
    required = ["filename", "calcmass", "rt", "charge"]
    missing = [c for c in required if conf._optional_columns[c] is None]
    if missing:
        missing = ", ".join([c + "_column" for c in missing])
        raise ValueError(
            "The following parameters must be specified when loading a "
            "collection of PSMs in order to save them in FlashLFQ format: "
            f"{missing}"
        )

    if conf._has_proteins:
        proteins = conf._proteins
    elif conf._protein_column is not None:
        proteins = conf._protein_column
    else:
        proteins = None

    # Get parameters
    peptides = conf.peptides
    filename_column = conf._optional_columns["filename"]
    peptide_column = conf._peptide_column
    mass_column = conf._optional_columns["calcmass"]
    rt_column = conf._optional_columns["rt"]
    charge_column = conf._optional_columns["charge"]
    eval_fdr = conf._eval_fdr

    # Create FlashLFQ dataframe
    logger.info("FDR threshold for FlashLFQ export: %.4f", eval_fdr)
    passing = peptides["mokapot q-value"] <= eval_fdr

    out_df = pd.DataFrame()
    out_df["File Name"] = peptides.loc[passing, filename_column].apply(lambda x: Path(x).name)

    seq = peptides.loc[passing, peptide_column]
    base_seq = (
        seq.str.replace(r"[\[\(].*?[\]\)]", "", regex=True)
        .str.replace(r"^.*?\.", "", regex=True)
        .str.replace(r"\..*?$", "", regex=True)
    )

    out_df["Base Sequence"] = base_seq
    out_df["Full Sequence"] = seq
    out_df["Peptide Monoisotopic Mass"] = peptides.loc[passing, mass_column]
    out_df["Scan Retention Time"] = peptides.loc[passing, rt_column] / 60
    out_df["Precursor Charge"] = peptides.loc[passing, charge_column]

    if isinstance(proteins, str):
        # TODO: Add delimiter sniffing.
        prots = peptides.loc[passing, proteins].str.replace("\t", "; ", regex=False)
    elif proteins is None:
        prots = ""
    else:
        prots = base_seq.map(proteins.peptide_map.get)
        shared = pd.isna(prots)
        prots.loc[shared] = base_seq[shared].map(proteins.shared_peptides.get)

    out_df["Protein Accession"] = prots
    missing = pd.isna(out_df["Protein Accession"])
    num_missing = missing.sum()
    if num_missing:
        logger.warning(
            "- Discarding %i peptides that could not be mapped to protein groups",
            num_missing,
        )
        out_df = out_df.loc[~missing, :]

    return out_df
