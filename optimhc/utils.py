# utils.py

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Optional, Dict, Union
from logging import getLogger

logger = getLogger(__name__)


def convert_pfm_to_pwm(pfm_filename, pseudocount=0.8, background_freqs=None):
    """
    Convert a Position Frequency Matrix (PFM) file to a Position Weight Matrix (PWM).

    Parameters
    ----------
    pfm_filename : str
        The file path to the PFM file.
    pseudocount : float, optional
        The pseudocount to add to the PFM to avoid zero probabilities. Default is 0.8.
    background_freqs : dict, optional
        Dictionary containing the background frequencies for each amino acid.
        If None, uses 1/20 for all.

    Returns
    -------
    pd.DataFrame
        DataFrame representation of the PWM.

    Notes
    -----
    The conversion process involves:
    1. Adding pseudocounts to the PFM
    2. Converting to Position Probability Matrix (PPM)
    3. Converting to PWM using log2(PPM/background_freqs)
    """
    # Default background frequencies if not provided
    if background_freqs is None:
        background_freqs = 1 / 20

    pfm = pd.read_csv(pfm_filename, sep="\t", header=None, index_col=0)
    pfm.drop(pfm.columns[-1], axis=1, inplace=True)  # Drop any extraneous columns
    pfm_pseudo = pfm + pseudocount
    ppm = pfm_pseudo.div(pfm_pseudo.sum(axis=0), axis=1)
    pwm = np.log2(ppm.div(list(background_freqs.values())))

    return pwm


def strip_flanking_and_charge(peptide: str) -> str:
    """
    Remove the pre and next amino acids from a peptide sequence.
    Also when there is a charge state at the end of the peptide, remove it.

    Parameters
    ----------
    peptide : str
        The peptide sequence with flanking amino acids.
        Example: '.AANDAGYFNDEMAPIEVKTK.'
        Example: 'F.VTVQGRAIC[119.0041]SDPNNKRVKN4.A'
        Example: '-.RRVEHHDHAVVSGR4.L'

    Returns
    -------
    str
        The peptide sequence with flanking amino acids removed.
        Example: 'AANDAGYFNDEMAPIEVKTK'
        Example: 'VTVQGRAIC[119.0041]SDPNNKRVKN'
        Example: 'RRVEHHDHAVVSGR'

    Notes
    -----
    This function removes any amino acids before the first '.' and after the last '.'
    in the peptide sequence.
    """
    import re

    peptide = re.sub(r"^[^.]*\.|\.[^.]*$", "", peptide)
    
    # Some PIN may have charge state at the end of the peptide, e.g., R.RRVEHHDHAVVSGR4.L
    # We should remove it to get the correct peptide sequence
    # Fragpipe: F.VTVQGRAIC[119.0041]SDPNNKRVKN4.A, remove charge state '4' in the peptide
    
    if peptide and peptide[-1].isdigit():
        peptide = peptide[:-1]
    
    return peptide
    

    


def remove_modifications(peptide: str, keep_modification=None) -> str:
    """
    Remove modifications from a peptide sequence, with an option to keep specific modifications.

    Parameters
    ----------
    peptide : str
        The peptide sequence with modifications in brackets.
        Example: 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    keep_modification : str or list, optional
        The modification(s) to keep. If provided, only these modifications will be
        preserved in the output sequence. Default is None.

    Returns
    -------
    str
        The peptide sequence with modifications removed or kept.
        Example: 'AANDAGYFNDEMAPIEVKTK' (if keep_modification is None)
        Example: 'AANDAGYFNDEM[15.9949]APIEVKTK' (if keep_modification=['15.9949'])

    Notes
    -----
    Modifications are specified in square brackets after the amino acid.
    If keep_modification is provided, only those specific modifications will be
    preserved in the output sequence.
    """
    import re

    if keep_modification is None:
        return re.sub(r"\[.*?\]", "", peptide)
    else:
        if not isinstance(keep_modification, list):
            keep_modification = [keep_modification]

        def replace_mod(match):
            mod = match.group(0)
            if any(keep in mod for keep in keep_modification):
                return mod
            return ""

        return re.sub(r"\[.*?\]", replace_mod, peptide)


def preprocess_peptide(peptide: str) -> str:
    """
    Preprocess the peptide sequence by removing flanking regions and modifications.

    Parameters
    ----------
    peptide : str
        Original peptide sequence with possible flanking regions and modifications.
        Example: '.AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK.'

    Returns
    -------
    str
        Cleaned peptide sequence without flanking regions and modifications.
        Example: 'AANDAGYFNDEMAPIEVKTK'

    Notes
    -----
    This function performs two operations in sequence:
    1. Removes flanking amino acids using remove_pre_and_nxt_aa
    2. Removes all modifications using remove_modifications
    """
    peptide = strip_flanking_and_charge(peptide)
    peptide = remove_modifications(peptide)
    return peptide


def list_all_files_in_directory(directory_path: str) -> List[str]:
    """
    Retrieve all files in the specified directory and return a list of file paths.

    Parameters
    ----------
    directory_path : str
        The path to the directory to search in.
        Example: '/path/to/directory'

    Returns
    -------
    list of str
        List of absolute file paths found in the directory and its subdirectories.
        Example: ['/path/to/directory/file1.txt', '/path/to/directory/subdir/file2.txt']

    Notes
    -----
    This function recursively searches through all subdirectories and returns
    absolute paths for all files found.
    """
    path = Path(directory_path)
    file_list = [str(file) for file in path.rglob("*") if file.is_file()]
    return file_list


def extract_unimod_from_peptidoform(peptide: str, mod_dict: dict) -> tuple:
    """
    Convert a modified peptide sequence into DeepLC format.

    Parameters
    ----------
    peptide : str
        The input peptide sequence with modifications in brackets.
        Example: 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    mod_dict : dict
        Dictionary mapping modification names (in peptide) to corresponding Unimod names.
        Example: {'15.9949': 'Oxidation', '42.0106': 'Acetyl'}

    Returns
    -------
    tuple
        (seq, modifications):
            seq : str
                The unmodified peptide sequence.
            modifications : str
                String of modifications formatted as `position|UnimodName`, separated by pipes `|`.
    """
    clean_sequence = ""
    modifications = []
    current_position = 0
    i = 0
    while i < len(peptide):
        if peptide[i] == "[":
            end = peptide.find("]", i)
            if end == -1:
                raise ValueError(
                    f"Invalid modification format in {peptide}: missing closing bracket."
                )
            mod_name = peptide[i + 1 : end]
            if mod_name not in mod_dict:
                raise ValueError(
                    f"Modification '{mod_name}' not found in the dictionary."
                )
            modifications.append(f"{current_position}|{mod_dict[mod_name]}")
            i = end + 1
        else:
            clean_sequence += peptide[i]
            current_position += 1
            i += 1

    modification_str = "|".join(modifications)
    logger.debug(f"Original peptide: {peptide}")
    logger.debug(f"Output clean_sequence: {clean_sequence}")
    logger.debug(f"Output modifications: {modification_str}")
    return clean_sequence, modification_str


def convert_to_unimod_format(peptide: str, mod_dict: dict) -> str:
    """
    Convert a modified peptide sequence into Unimod format.

    Parameters
    ----------
    peptide : str
        The input peptide sequence with modifications in brackets.
        Example: 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    mod_dict : dict
        Dictionary mapping modification names (in peptide) to corresponding Unimod names.
        Example: {'15.9949': 'UNIMOD:4', '42.0106': 'UNIMOD:1'}

    Returns
    -------
    str
        The peptide sequence formatted for Unimod.
        Example: 'AANDAGYFNDEM[UNIMOD:4]APIEVK[UNIMOD:1]TK'

    Notes
    -----
    This function replaces the modification names in brackets with their
    corresponding Unimod identifiers while preserving the peptide sequence
    structure.
    """
    res = peptide
    for key, value in mod_dict.items():
        res = res.replace(f"[{key}]", f"[{value}]")
    logger.debug(f"Original peptide: {peptide}")
    logger.debug(f"Output peptide: {res}")
    return res
