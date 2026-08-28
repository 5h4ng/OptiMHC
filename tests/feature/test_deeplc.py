import importlib
import sys
from types import ModuleType
from unittest.mock import patch

import numpy as np
import pandas as pd

from optimhc.psm_container import PsmContainer


class RecordingDeepLC:
    calibration_calls = []
    prediction_calls = []

    def __init__(self, **kwargs):
        self.calibrated_rt = None

    def calibrate_preds(self, *, seq_df):
        self.calibration_calls.append(seq_df.copy())
        self.calibrated_rt = float(seq_df["tr"].iloc[0])

    def make_preds(self, *, seq_df):
        self.prediction_calls.append(seq_df.copy())
        return np.full(len(seq_df), self.calibrated_rt)


deeplc_stub = ModuleType("deeplc")
deeplc_stub.DeepLC = RecordingDeepLC
with patch.dict(sys.modules, {"deeplc": deeplc_stub}):
    deeplc_module = importlib.import_module("optimhc.feature.deeplc")
DeepLCFeatureGenerator = deeplc_module.DeepLCFeatureGenerator


def test_deeplc_calibrates_and_predicts_each_acquisition_run_separately():
    RecordingDeepLC.calibration_calls = []
    RecordingDeepLC.prediction_calls = []
    psms = PsmContainer(
        pd.DataFrame(
            {
                "psm_id": [0, 1, 2, 3],
                "run": ["run_a", "run_b", "run_a", "run_b"],
                "scan": [10, 10, 11, 11],
                "rank": [1, 1, 1, 1],
                "sequence": ["AAAA", "BBBB", "CCCC", "DDDD"],
                "mods": ["", "", "", ""],
                "mod_sites": ["", "", "", ""],
                "charge": [2, 2, 2, 2],
                "proteins": ["P1", "P2", "P3", "P4"],
                "is_decoy": [False, False, False, False],
                "retention_time": [10.0, 100.0, 20.0, 200.0],
                "search_score": [4.0, 3.0, 2.0, 1.0],
            }
        ),
        feature_columns=("search_score",),
    )
    generator = DeepLCFeatureGenerator(
        psms,
        calibration_criteria_column="search_score",
        calibration_set_size=1,
    )

    features = generator.generate_features()

    assert [call["seq"].tolist() for call in RecordingDeepLC.calibration_calls] == [
        ["AAAA"],
        ["BBBB"],
    ]
    assert [call["seq"].tolist() for call in RecordingDeepLC.prediction_calls] == [
        ["AAAA", "CCCC"],
        ["BBBB", "DDDD"],
    ]
    assert features["predicted_retention_time"].tolist() == [10.0, 100.0, 10.0, 100.0]
