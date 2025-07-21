import pytest
import pandas as pd
import numpy as np
import tempfile
import os

from optimhc.psm_container import PsmContainer


class TestPsmContainerBasics:
    """Basic functionality tests for PsmContainer - these tests should not change in future iterations"""
    
    def test_initialization_success(self, minimal_psm_data, basic_rescoring_features):
        """Test successful initialization"""
        container = PsmContainer(
            psms=minimal_psm_data,
            label_column='label',
            scan_column='scan_id',
            spectrum_column='spectrum_id',
            ms_data_file_column='ms_file',
            peptide_column='peptide',
            protein_column='protein',
            rescoring_features=basic_rescoring_features
        )
        
        assert len(container) == 4
        assert container.label_column == 'label'
        assert container.scan_column == 'scan_id'
        assert container.spectrum_column == 'spectrum_id'
        assert container.ms_data_file_column == 'ms_file'
        assert container.peptide_column == 'peptide'
        assert container.protein_column == 'protein'
        assert container.rescoring_features == basic_rescoring_features

    def test_initialization_missing_column(self, minimal_psm_data, basic_rescoring_features):
        """Test initialization failure when required column is missing"""
        with pytest.raises(ValueError, match="Column 'missing_column' not found"):
            PsmContainer(
                psms=minimal_psm_data,
                label_column='missing_column',
                scan_column='scan_id',
                spectrum_column='spectrum_id',
                ms_data_file_column='ms_file',
                peptide_column='peptide',
                protein_column='protein',
                rescoring_features=basic_rescoring_features
            )

    def test_initialization_no_decoys(self, basic_rescoring_features):
        """Test initialization failure when no decoy PSMs are present"""
        all_target_data = pd.DataFrame({
            'scan_id': [1, 2, 3],
            'spectrum_id': ['spec_1', 'spec_2', 'spec_3'],
            'peptide': ['PEPTIDE1', 'PEPTIDE2', 'PEPTIDE3'],
            'protein': ['PROT1', 'PROT2', 'PROT3'],
            'label': [True, True, True],
            'ms_file': ['file1.mzML', 'file1.mzML', 'file2.mzML'],
            'feature1': [0.1, 0.2, 0.3],
            'feature2': [1.1, 1.2, 1.3]
        })
        
        with pytest.raises(ValueError, match="No decoy PSMs found"):
            PsmContainer(
                psms=all_target_data,
                label_column='label',
                scan_column='scan_id',
                spectrum_column='spectrum_id',
                ms_data_file_column='ms_file',
                peptide_column='peptide',
                protein_column='protein',
                rescoring_features=basic_rescoring_features
            )

    def test_initialization_no_targets(self, basic_rescoring_features):
        """Test initialization failure when no target PSMs are present"""
        all_decoy_data = pd.DataFrame({
            'scan_id': [1, 2, 3],
            'spectrum_id': ['spec_1', 'spec_2', 'spec_3'],
            'peptide': ['DECOY_PEPTIDE1', 'DECOY_PEPTIDE2', 'DECOY_PEPTIDE3'],
            'protein': ['DECOY_PROT1', 'DECOY_PROT2', 'DECOY_PROT3'],
            'label': [False, False, False],
            'ms_file': ['file1.mzML', 'file1.mzML', 'file2.mzML'],
            'feature1': [0.1, 0.2, 0.3],
            'feature2': [1.1, 1.2, 1.3]
        })
        
        with pytest.raises(ValueError, match="All PSMs are labeled as decoy"):
            PsmContainer(
                psms=all_decoy_data,
                label_column='label',
                scan_column='scan_id',
                spectrum_column='spectrum_id',
                ms_data_file_column='ms_file',
                peptide_column='peptide',
                protein_column='protein',
                rescoring_features=basic_rescoring_features
            )

    def test_basic_properties(self, psm_container):
        """Test basic property access"""
        assert len(psm_container) == 4
        assert len(psm_container.target_psms) == 3
        assert len(psm_container.decoy_psms) == 1
        assert len(psm_container.peptides) == 4
        assert set(psm_container.feature_columns) == {'feature1', 'feature2'}
        assert psm_container.feature_sources == ['basic']

    def test_psms_property_returns_copy(self, psm_container):
        """Test that psms property returns a copy, not affecting original data"""
        original_len = len(psm_container)
        psms_copy = psm_container.psms
        
        # Modifying the copy should not affect the original data
        psms_copy.loc[0, 'peptide'] = 'MODIFIED'
        
        assert len(psm_container) == original_len
        assert psm_container.peptides[0] != 'MODIFIED'

    def test_target_decoy_separation(self, psm_container):
        """Test correct separation of target and decoy PSMs"""
        targets = psm_container.target_psms
        decoys = psm_container.decoy_psms
        
        # Check target PSMs
        assert all(targets['label'] == True)
        assert len(targets) == 3
        
        # Check decoy PSMs
        assert all(decoys['label'] == False)
        assert len(decoys) == 1
        
        # Check total count
        assert len(targets) + len(decoys) == len(psm_container)

    def test_copy_method(self, psm_container):
        """Test deep copy method"""
        copied_container = psm_container.copy()
        
        # Check they are different objects
        assert copied_container is not psm_container
        assert copied_container._psms is not psm_container._psms
        
        # Check content is the same
        assert len(copied_container) == len(psm_container)
        assert copied_container.rescoring_features == psm_container.rescoring_features
        pd.testing.assert_frame_equal(copied_container.psms, psm_container.psms)

    def test_repr_method(self, psm_container):
        """Test string representation method"""
        repr_str = repr(psm_container)
        
        assert "PsmContainer with 4 PSMs" in repr_str
        assert "Target PSMs: 3" in repr_str
        assert "Decoy PSMs: 1" in repr_str
        assert "Unique Peptides:" in repr_str
        assert "Rescoring Features:" in repr_str

    def test_len_method(self, psm_container):
        """Test length method"""
        assert len(psm_container) == 4


class TestPsmContainerFeatureManagement:
    """Tests for feature management functionality"""
    
    def test_add_features_basic(self, psm_container):
        """Test basic feature addition"""
        new_features = pd.DataFrame({
            'scan_id': [1, 2, 3, 4],
            'new_feature1': [0.5, 0.6, 0.7, 0.8],
            'new_feature2': [2.1, 2.2, 2.3, 2.4]
        })
        
        original_len = len(psm_container)
        psm_container.add_features(
            features_df=new_features,
            psms_key='scan_id',
            feature_key='scan_id',
            source='new_source'
        )
        
        # Check length remains unchanged
        assert len(psm_container) == original_len
        
        # Check new features are added
        assert 'new_feature1' in psm_container.columns
        assert 'new_feature2' in psm_container.columns
        
        # Check feature source is correctly recorded
        assert 'new_source' in psm_container.rescoring_features
        assert set(psm_container.rescoring_features['new_source']) == {'new_feature1', 'new_feature2'}

    def test_add_features_with_suffix(self, psm_container):
        """Test feature addition with suffix"""
        new_features = pd.DataFrame({
            'scan_id': [1, 2, 3, 4],
            'feature1': [0.9, 1.0, 1.1, 1.2]  # Same name as existing feature
        })
        
        psm_container.add_features(
            features_df=new_features,
            psms_key='scan_id',
            feature_key='scan_id',
            source='conflicting_source',
            suffix='_new'
        )
        
        # Check suffixed feature is added
        assert 'feature1_new' in psm_container.columns
        assert 'feature1' in psm_container.columns  # Original feature still exists
        
        # Check feature source record
        assert 'feature1_new' in psm_container.rescoring_features['conflicting_source']

    def test_drop_features(self, psm_container):
        """Test feature deletion"""
        original_features = psm_container.feature_columns.copy()
        
        psm_container.drop_features(['feature1'])
        
        # Check feature is deleted
        assert 'feature1' not in psm_container.columns
        assert 'feature2' in psm_container.columns
        
        # Check feature source is updated
        assert 'feature1' not in psm_container.rescoring_features['basic']
        assert 'feature2' in psm_container.rescoring_features['basic']

    def test_drop_source(self, psm_container):
        """Test feature source deletion"""
        psm_container.drop_source('basic')
        
        # Check all related features are deleted
        assert 'feature1' not in psm_container.columns
        assert 'feature2' not in psm_container.columns
        
        # Check feature source is deleted
        assert 'basic' not in psm_container.rescoring_features


class TestPsmContainerDataExport:
    """Tests for data export functionality"""
    
    def test_write_pin_basic(self, psm_container):
        """Test basic PIN file writing functionality"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pin', delete=False) as tmp_file:
            tmp_path = tmp_file.name
        
        try:
            result_df = psm_container.write_pin(tmp_path)
            
            # Check file is created
            assert os.path.exists(tmp_path)
            
            # Check DataFrame structure
            expected_columns = ['SpecId', 'Label', 'ScanNr', 'feature1', 'feature2', 'Peptide', 'Proteins']
            assert list(result_df.columns) == expected_columns
            
            # Check label conversion
            assert set(result_df['Label'].unique()) == {1, -1}
            assert len(result_df) == len(psm_container)
            
            # Read file to verify format
            saved_df = pd.read_csv(tmp_path, sep='\t')
            pd.testing.assert_frame_equal(result_df, saved_df)
            
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)

    def test_write_pin_with_specific_source(self, psm_container):
        """Test PIN file writing with specific feature source"""
        # Add new feature source
        new_features = pd.DataFrame({
            'scan_id': [1, 2, 3, 4],
            'extra_feature': [0.5, 0.6, 0.7, 0.8]
        })
        
        psm_container.add_features(
            features_df=new_features,
            psms_key='scan_id',
            feature_key='scan_id',
            source='extra'
        )
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pin', delete=False) as tmp_file:
            tmp_path = tmp_file.name
        
        try:
            # Only export features from basic source
            result_df = psm_container.write_pin(tmp_path, source=['basic'])
            
            # Check only specified source features are included
            assert 'feature1' in result_df.columns
            assert 'feature2' in result_df.columns
            assert 'extra_feature' not in result_df.columns
            
        finally:
            if os.path.exists(tmp_path):
                os.unlink(tmp_path)


class TestPsmContainerAdvancedFeatures:
    """Tests for advanced PsmContainer features"""
    
    def test_get_top_hits_with_rank(self, complex_psm_container):
        """Test getting top hits based on rank"""
        top1 = complex_psm_container.get_top_hits(n=1)
        
        # top 1 should only include PSMs with hit_rank <= 1
        assert all(top1.psms['hit_rank'] <= 1)
        
        # Check that we have the expected number of top hits
        expected_top1_count = sum(complex_psm_container.psms['hit_rank'] <= 1)
        assert len(top1) == expected_top1_count

    def test_get_top_hits_without_rank(self, minimal_psm_data, basic_rescoring_features):
        """Test get_top_hits when rank column is not specified"""
        container = PsmContainer(
            psms=minimal_psm_data,
            label_column='label',
            scan_column='scan_id',
            spectrum_column='spectrum_id',
            ms_data_file_column='ms_file',
            peptide_column='peptide',
            protein_column='protein',
            rescoring_features=basic_rescoring_features
            # No hit_rank_column specified
        )
        
        top1 = container.get_top_hits(n=1)
        assert len(top1) == len(container)  # Should return all PSMs


if __name__ == "__main__":
    pytest.main([__file__]) 