import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from optimhc.parser import read_pepxml, read_pin

DEFAULT_DATA = Path(__file__).parents[2] / "data" / "psm-container-refactor"
REAL_DATA = Path(os.environ.get("OPTIMHC_TEST_DATA", DEFAULT_DATA))
pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not REAL_DATA.exists(), reason="local real-data fixtures are unavailable"),
]


def test_msfragger_pin_and_pepxml_have_the_same_candidate_identity():
    stem = "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn01"
    pin = read_pin(REAL_DATA / f"{stem}.pin")
    pepxml = read_pepxml(REAL_DATA / f"{stem}.pepXML", decoy_prefix="rev_")

    identity = [
        "run",
        "scan",
        "rank",
        "sequence",
        "mods",
        "mod_sites",
        "charge",
        "is_decoy",
    ]
    pd.testing.assert_frame_equal(pin.df[identity], pepxml.df[identity])
    np.testing.assert_allclose(
        pin.df["retention_time"],
        pepxml.df["retention_time"],
        atol=1e-3,
    )


def test_msbooster_pin_files_keep_their_runs_and_features():
    containers = [
        read_pin(REAL_DATA / f"LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn0{run}_edited.pin")
        for run in (1, 2)
    ]

    assert [len(container) for container in containers] == [10_515, 20_017]
    assert [container.df["run"].nunique() for container in containers] == [1, 1]
    assert containers[0].df["run"].iat[0] != containers[1].df["run"].iat[0]
    assert {
        "pred_RT_real_units",
        "unweighted_spectral_entropy",
        "delta_RT_loess",
    }.issubset(containers[0].feature_columns)


def test_pipeline_combines_msbooster_runs_before_assigning_psm_ids(tmp_path):
    from optimhc.core.pipeline import Pipeline

    inputs = [
        str(REAL_DATA / f"LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn0{run}_edited.pin")
        for run in (1, 2)
    ]
    pipeline = Pipeline(
        {
            "inputType": "pin",
            "inputFile": inputs,
            "outputDir": str(tmp_path),
            "visualization": False,
            "saveModels": False,
        }
    )

    candidates = pipeline.read_input()

    assert len(candidates) == 30_532
    assert candidates.df["psm_id"].tolist() == list(range(30_532))
    assert candidates.df["run"].drop_duplicates().tolist() == [
        "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn01",
        "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn02",
    ]


def test_comet_pepxml_preserves_multi_rank_candidates():
    stem = "20200804_COL011-00022-HLAIp_R03"
    candidates = read_pepxml(REAL_DATA / f"{stem}.pep.xml", decoy_prefix="DECOY_")

    assert len(candidates) == 152_201
    assert set(candidates.df["rank"]) == {1, 2, 3, 4, 5}
    assert candidates.df["run"].unique().tolist() == [stem]
    assert (REAL_DATA / f"{stem}.mzML").exists()


def test_pipeline_combines_two_pepxml_runs_in_configured_order(tmp_path):
    from optimhc.core.pipeline import Pipeline

    inputs = [
        str(REAL_DATA / f"LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn0{run}.pepXML")
        for run in (1, 2)
    ]
    pipeline = Pipeline(
        {
            "inputType": "pepxml",
            "inputFile": inputs,
            "outputDir": str(tmp_path),
            "decoyPrefix": "rev_",
            "visualization": False,
            "saveModels": False,
        }
    )

    candidates = pipeline.read_input()

    assert len(candidates) == 30_532
    assert candidates.df["psm_id"].tolist() == list(range(30_532))
    assert candidates.df["run"].drop_duplicates().tolist() == [
        "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn01",
        "LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn02",
    ]


@pytest.mark.raw_data
@pytest.mark.skipif(
    os.environ.get("OPTIMHC_RUN_RAW_DATA") != "1",
    reason="set OPTIMHC_RUN_RAW_DATA=1 to parse the large mzML fixtures",
)
def test_two_run_mzml_lookup_joins_by_run_and_scan_without_spec_idx():
    from optimhc.parser import extract_mzml_data

    candidate_rows = []
    spectra = []
    for run_number in (1, 2):
        stem = f"LP20210421_HT_DCfeeding_Exp1_HLAII_0hr_fxn0{run_number}"
        candidate = read_pepxml(REAL_DATA / f"{stem}.pepXML", decoy_prefix="rev_").df.iloc[0]
        candidate_rows.append({"run": stem, "scan": candidate["scan"]})
        spectra.append(
            extract_mzml_data(
                str(REAL_DATA / f"{stem}.mzML"),
                scan_ids=[int(candidate["scan"])],
            ).rename(columns={"source": "run"})
        )

    joined = pd.DataFrame(candidate_rows).merge(
        pd.concat(spectra, ignore_index=True),
        on=["run", "scan"],
        validate="one_to_one",
    )

    assert len(joined) == 2
    assert joined["mz"].map(len).gt(0).all()
    assert joined["intensity"].map(len).gt(0).all()
