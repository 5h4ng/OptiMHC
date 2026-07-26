"""Tests for spectral alignment: exp/pred spectra -> two aligned intensity vectors."""

import numpy as np
import pandas as pd
import pytest

from optimhc.feature.spectral_similarity import SpectralSimilarityFeatureGenerator


def _make_generator(tolerance_ppm=20.0, top_n=36):
    gen = object.__new__(SpectralSimilarityFeatureGenerator)
    gen.tolerance_ppm = tolerance_ppm
    gen.top_n = top_n
    return gen


PRED_MZ = [147.11, 234.15, 347.23, 462.26, 577.29, 690.37, 803.45, 916.54]
PRED_INT = [0.05, 0.30, 0.80, 1.00, 0.60, 0.35, 0.15, 0.02]

EXP_MZ = [
    130.07,
    147.112,
    234.153,
    300.10,
    347.234,
    462.265,
    520.33,
    577.30,
    690.38,
    750.40,
    803.46,
    950.00,
]
EXP_INT = [
    500.0,
    1200.0,
    8000.0,
    300.0,
    15000.0,
    22000.0,
    400.0,
    12000.0,
    7500.0,
    200.0,
    3000.0,
    100.0,
]


class TestAlignSpectraAllPeaks:
    def test_realistic_alignment(self):
        gen = _make_generator(tolerance_ppm=20.0)
        aligned_exp, aligned_pred, matched, info = gen._align_spectra_all_peaks(
            EXP_MZ,
            EXP_INT,
            PRED_MZ,
            PRED_INT,
        )
        assert len(aligned_exp) == len(PRED_MZ)
        assert len(aligned_pred) == len(PRED_MZ)
        assert matched.shape == (len(PRED_MZ), 2)

        # pred 147.11 should match exp 147.112
        assert aligned_exp[0] == 1200.0
        # pred 462.26 should match exp 462.265
        assert aligned_exp[3] == 22000.0
        # pred 916.54 has no nearby exp peak -> 0
        assert aligned_exp[7] == 0.0
        # unmatched index should be -1
        assert matched[7, 1] == -1


class TestGetTopPeaksVectors:
    def test_top_n_from_aligned(self):
        gen = _make_generator(tolerance_ppm=20.0)
        aligned_exp, aligned_pred, matched, _ = gen._align_spectra_all_peaks(
            EXP_MZ,
            EXP_INT,
            PRED_MZ,
            PRED_INT,
        )
        top_exp, top_pred, top_matched = gen._get_top_peaks_vectors(
            aligned_exp,
            aligned_pred,
            matched,
            top_n=3,
        )
        assert len(top_exp) == 3
        assert top_matched.shape == (3, 2)
        # top 3 by pred intensity: 462.26 (1.0), 347.23 (0.8), 577.29 (0.6)
        np.testing.assert_array_equal(sorted(top_pred), sorted([1.0, 0.8, 0.6]))


class TestAlignPipeline:
    def test_e2e_align_and_top_peaks(self):
        gen = _make_generator(tolerance_ppm=20.0, top_n=5)
        aligned_exp, aligned_pred, matched, _ = gen._align_spectra_all_peaks(
            EXP_MZ,
            EXP_INT,
            PRED_MZ,
            PRED_INT,
        )
        top_exp, top_pred, _ = gen._get_top_peaks_vectors(
            aligned_exp,
            aligned_pred,
            matched,
            top_n=5,
        )
        features = gen._calculate_similarity_features(top_exp, top_pred)
        assert features["spectral_angle_similarity"] > 0.5
        assert features["predicted_seen_nonzero"] >= 4


def test_generate_features_requires_experimental_spectra(monkeypatch):
    gen = _make_generator()
    gen.df = pd.DataFrame(
        {
            "processed_peptide": ["PEPTIDE"],
            "charge": [2],
            "scan": [10],
            "mz_file_path": ["run.mzML"],
        }
    )
    monkeypatch.setattr(
        gen,
        "_predict_theoretical_spectra",
        lambda processed_peptides, charges: pd.DataFrame(),
    )
    monkeypatch.setattr(gen, "_extract_experimental_spectra", pd.DataFrame)

    with pytest.raises(ValueError, match="could not extract experimental spectra"):
        gen._generate_features()
