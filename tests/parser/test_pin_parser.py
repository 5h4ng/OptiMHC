import pytest
import pandas as pd
import os
from pathlib import Path
import re

from optimhc.parser.pin import read_pin, _read_single_pin_as_df
from optimhc.psm_container import PsmContainer


class TestPinParser:
    """Test the PIN parser functionality"""

    @pytest.fixture
    def fragpipe_pin_path(self):
        """Return the path to the test FragPipe PIN file"""
        return os.path.join(os.path.dirname(__file__), "fragpipe_sample.pin")
    
    @pytest.fixture
    def percolator_pin_path(self):
        """Return the path to the test Percolator PIN file"""
        return os.path.join(os.path.dirname(__file__), "percolator_sample.pin")

    def test_read_single_pin_fragpipe(self, fragpipe_pin_path):
        """Test reading a single FragPipe PIN file"""
        df = _read_single_pin_as_df(fragpipe_pin_path)
        
        # Check that the DataFrame has the correct structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 4  # Should have 4 PSMs
        
        expected_columns = ["SpecId", "Label", "ScanNr", "ExpMass", "retentiontime", "rank", 
                           "Peptide", "Proteins"]
        for col in expected_columns:
            assert col in df.columns
        
        # Check target PSMs
        assert df.iloc[0]["SpecId"] == "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn01.5401.5401.4_1"
        assert df.iloc[0]["Label"] == "1"
        assert df.iloc[0]["ScanNr"] == "5401"
        assert df.iloc[0]["Peptide"] == "R.RRVEHHDHAVVSGR4.L"
        
        # Check decoy PSMs
        assert df.iloc[3]["Label"] == "-1"
        assert "rev_" in df.iloc[3]["Proteins"] or df.iloc[3]["Proteins"].startswith("rev_")

    def test_read_single_pin_percolator(self, percolator_pin_path):
        """Test reading a single Percolator PIN file"""
        df = _read_single_pin_as_df(percolator_pin_path)
        
        expected_columns = ["SpecId", "Label", "ScanNr", "ExpMass", "CalcMass", 
                           "lnExpect", "Xcorr", "Peptide", "Proteins"]
        for col in expected_columns:
            assert col in df.columns
        
        assert df.iloc[0]["Label"] == "-1"  
        assert "ZBTB4_HUMAN" in df.iloc[0]["Proteins"]
        assert df.iloc[1]["Label"] == "1"  

    def test_read_pin_fragpipe_to_container(self, fragpipe_pin_path):
        """Test reading a FragPipe PIN file into a PsmContainer"""
        container = read_pin(fragpipe_pin_path, retention_time_column="retentiontime")
        
        # Check that a PsmContainer was created
        assert isinstance(container, PsmContainer)
        
        # Check that the PsmContainer has the correct data
        assert len(container.psms) == 4
        
        # Check charge extraction from charge_* columns
        assert "Charge" in container.psms.columns
        # In the sample data, the first PSM has charge_4=1, so charge should be 4
        assert container.psms.iloc[0]["Charge"] == 4
        assert container.psms.iloc[3]["Charge"] == 5  # Last PSM has charge_5=1
        
        # Check that the target/decoy labels were properly converted to boolean
        assert container.psms[container.label_column].dtype == bool
        assert container.psms.iloc[0][container.label_column] == True   
        assert container.psms.iloc[1][container.label_column] == True   
        assert container.psms.iloc[2][container.label_column] == True   
        assert container.psms.iloc[3][container.label_column] == False  
        
        # Check that spectrum IDs were parsed to extract the hit rank
        assert container.psms[container.hit_rank_column].dtype == int
        assert all(container.psms[container.hit_rank_column] == 1)  # All rank 1 in the sample
        
        # Check retention time
        assert container.retention_time_column == "retentiontime"
        assert container.psms[container.retention_time_column].dtype == float

    def test_read_pin_percolator_to_container(self, percolator_pin_path):
        """Test reading a Percolator PIN file into a PsmContainer"""
        container = read_pin(percolator_pin_path)
        
        # Check that a PsmContainer was created
        assert isinstance(container, PsmContainer)
        
        # Check charge extraction from charge_* columns
        assert "Charge" in container.psms.columns
        # In percolator, the charge should be properly extracted from Charge1, Charge2, etc.
        assert container.psms.iloc[0]["Charge"] == 2
        assert container.psms.iloc[1]["Charge"] == 2
        assert container.psms.iloc[2]["Charge"] == 3
        
        # Check that the target/decoy labels were properly converted to boolean
        assert container.psms[container.label_column].dtype == bool
        assert not container.psms.iloc[0][container.label_column]  
        assert container.psms.iloc[1][container.label_column]      
        assert not container.psms.iloc[2][container.label_column]  
        
    def test_parse_specid(self, fragpipe_pin_path, percolator_pin_path):
        """Test the parsing of SpecId to extract hit ranks"""
        
        container_fp = read_pin(fragpipe_pin_path) 
        assert "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn01.5401.5401.4_1" in container_fp.psms["SpecId"].values
        
        container_perc = read_pin(percolator_pin_path)
        assert "/media/kevin/Elements/MhcValidator_analysis/data/JY_serial_dilution/IP0040_11MAI2022_JY_MHC1_HUMAN_S4_2PELLETS_KK_SERIAL_DIL_PT_3_R2_3_2_1" in container_perc.psms["SpecId"].values


