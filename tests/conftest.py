import tempfile
from pathlib import Path

import pandas as pd
import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files"""
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture
def sample_mzml_files(temp_dir):
    """Create mock mzML file paths"""
    files = []
    for i in range(3):
        file_path = Path(temp_dir) / f"sample_{i}.mzML"
        file_path.touch()
        files.append(str(file_path))
    return files


@pytest.fixture
def minimal_psm_data():
    """Create minimal PSM data for basic tests"""
    return pd.DataFrame(
        {
            "scan_id": [1, 2, 3, 4],
            "spectrum_id": ["spec_1", "spec_2", "spec_3", "spec_4"],
            "peptide": ["PEPTIDE1", "PEPTIDE2", "PEPTIDE3", "DECOY_PEPTIDE1"],
            "protein": ["PROT1", "PROT2", "PROT3", "DECOY_PROT1"],
            "label": [True, True, True, False],
            "ms_file": ["file1.mzML", "file1.mzML", "file2.mzML", "file1.mzML"],
            "feature1": [0.1, 0.2, 0.3, 0.4],
            "feature2": [1.1, 1.2, 1.3, 1.4],
        }
    )


@pytest.fixture
def complex_psm_data():
    """Create complex PSM data with additional columns for advanced tests"""
    return pd.DataFrame(
        {
            "scan_id": [1, 2, 3, 4, 5, 6],
            "spectrum_id": ["spec_1", "spec_2", "spec_3", "spec_4", "spec_5", "spec_6"],
            "peptide": [
                "PEPTIDE1",
                "PEPTIDE2",
                "PEPTIDE3",
                "DECOY_PEPTIDE1",
                "DECOY_PEPTIDE2",
                "PEPTIDE4",
            ],
            "protein": [
                "PROT1",
                "PROT2",
                "PROT3",
                "DECOY_PROT1",
                "DECOY_PROT2",
                "PROT4",
            ],
            "label": [True, True, True, False, False, True],
            "ms_file": [
                "file1.mzML",
                "file1.mzML",
                "file2.mzML",
                "file1.mzML",
                "file2.mzML",
                "file2.mzML",
            ],
            "charge": [2, 3, 2, 2, 3, 2],
            "hit_rank": [1, 1, 2, 1, 1, 1],
            "feature1": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            "feature2": [1.1, 1.2, 1.3, 1.4, 1.5, 1.6],
        }
    )


@pytest.fixture
def basic_rescoring_features():
    """Basic rescoring feature configuration"""
    return {"basic": ["feature1", "feature2"]}


@pytest.fixture
def complex_rescoring_features():
    """Complex rescoring feature configuration"""
    return {"basic": ["feature1", "feature2"], "advanced": ["feature3"]}


@pytest.fixture
def psm_container(minimal_psm_data, basic_rescoring_features):
    """Create a basic PsmContainer instance for testing"""
    from optimhc.psm_container import PsmContainer

    return PsmContainer(
        psms=minimal_psm_data,
        label_column="label",
        scan_column="scan_id",
        spectrum_column="spectrum_id",
        ms_data_file_column="ms_file",
        peptide_column="peptide",
        protein_column="protein",
        rescoring_features=basic_rescoring_features,
    )


@pytest.fixture
def complex_psm_container(complex_psm_data, basic_rescoring_features):
    """Create a complex PsmContainer instance with additional columns"""
    from optimhc.psm_container import PsmContainer

    return PsmContainer(
        psms=complex_psm_data,
        label_column="label",
        scan_column="scan_id",
        spectrum_column="spectrum_id",
        ms_data_file_column="ms_file",
        peptide_column="peptide",
        protein_column="protein",
        rescoring_features=basic_rescoring_features,
        hit_rank_column="hit_rank",
        charge_column="charge",
    )


def assert_dataframes_equal(df1, df2, check_dtype=True):
    """Helper function: Compare if two DataFrames are equal"""
    pd.testing.assert_frame_equal(df1, df2, check_dtype=check_dtype)


def assert_containers_equal(container1, container2):
    """Helper function: Compare if two PsmContainers are equal"""
    assert len(container1) == len(container2)
    assert container1.rescoring_features == container2.rescoring_features
    assert container1.label_column == container2.label_column
    assert container1.scan_column == container2.scan_column
    assert container1.spectrum_column == container2.spectrum_column
    pd.testing.assert_frame_equal(container1.psms, container2.psms)
