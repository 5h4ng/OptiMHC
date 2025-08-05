import pytest
import os
from optimhc import utils

def test_strip_flanking_and_charge():
    assert utils.strip_flanking_and_charge('.AANDAGYFNDEMAPIEVKTK.') == 'AANDAGYFNDEMAPIEVKTK'
    assert utils.strip_flanking_and_charge('F.VTVQGRAIC[119.0041]SDPNNKRVKN4.A') == 'VTVQGRAIC[119.0041]SDPNNKRVKN'
    assert utils.strip_flanking_and_charge('-.RRVEHHDHAVVSGR4.L') == 'RRVEHHDHAVVSGR'

def test_remove_modifications():
    pep = 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    assert utils.remove_modifications(pep) == 'AANDAGYFNDEMAPIEVKTK'
    assert utils.remove_modifications(pep, keep_modification='15.9949') == 'AANDAGYFNDEM[15.9949]APIEVKTK'
    assert utils.remove_modifications(pep, keep_modification=['42.0106']) == 'AANDAGYFNDEMAPIEVK[42.0106]TK'

def test_preprocess_peptide():
    pep = '.AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK.'
    assert utils.preprocess_peptide(pep) == 'AANDAGYFNDEMAPIEVKTK'

def test_list_all_files_in_directory(tmp_path):
    # Create files and subdirs
    d = tmp_path / "sub"
    d.mkdir()
    f1 = tmp_path / "file1.txt"
    f1.write_text("test")
    f2 = d / "file2.txt"
    f2.write_text("test2")
    files = utils.list_all_files_in_directory(str(tmp_path))
    assert str(f1) in files
    assert str(f2) in files

def test_extract_unimod_from_peptidoform():
    pep = 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    mod_dict = {'15.9949': 'Oxidation', '42.0106': 'Acetyl'}
    seq, mods = utils.extract_unimod_from_peptidoform(pep, mod_dict)
    assert seq == 'AANDAGYFNDEMAPIEVKTK'
    assert mods == '12|Oxidation|18|Acetyl'

def test_convert_to_unimod_format():
    pep = 'AANDAGYFNDEM[15.9949]APIEVK[42.0106]TK'
    mod_dict = {'15.9949': 'UNIMOD:4', '42.0106': 'UNIMOD:1'}
    res = utils.convert_to_unimod_format(pep, mod_dict)
    assert res == 'AANDAGYFNDEM[UNIMOD:4]APIEVK[UNIMOD:1]TK'
