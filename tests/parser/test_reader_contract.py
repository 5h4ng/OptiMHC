from pathlib import Path

import pandas as pd
import pytest

from optimhc.parser import read_pepxml, read_pin
from optimhc.rescore.mokapot import write_pin

FIXTURES = Path(__file__).parent


def test_pin_reader_returns_canonical_candidates():
    candidates = read_pin(FIXTURES / "fragpipe_sample.pin")

    assert candidates.df["psm_id"].tolist() == [0, 1, 2, 3]
    assert candidates.df["run"].nunique() == 1
    assert candidates.df["scan"].tolist() == [5401, 5459, 6277, 6587]
    assert candidates.df["rank"].tolist() == [1, 1, 1, 1]
    assert candidates.df["sequence"].tolist() == [
        "RRVEHHDHAVVSGR",
        "RQRPTRPPRQDKPPR",
        "RAAEDDEDDDVDTK",
        "QTDLMMVVWQEATLR",
    ]
    assert candidates.df["is_decoy"].tolist() == [False, False, False, True]
    assert candidates.df.loc[3, "mods"] == "Oxidation@M"
    assert candidates.df.loc[3, "mod_sites"] == "5"
    assert "rank" in candidates.feature_columns
    assert "SpecId" not in candidates.df.columns


def test_pepxml_reader_preserves_every_candidate_rank():
    candidates = read_pepxml(FIXTURES / "multi_rank_sample.pep.xml")

    assert candidates.df["psm_id"].tolist() == [0, 1]
    assert candidates.df["run"].tolist() == ["run_a", "run_a"]
    assert candidates.df["scan"].tolist() == [10, 10]
    assert candidates.df["rank"].tolist() == [1, 2]
    assert candidates.df["sequence"].tolist() == ["PEPM", "EDIT"]
    assert candidates.df["mods"].tolist() == ["Oxidation@M", ""]
    assert candidates.df["mod_sites"].tolist() == ["4", ""]
    assert candidates.df["proteins"].tolist() == ["P1", "DECOY_P1;DECOY_P2"]
    assert candidates.df["is_decoy"].tolist() == [False, True]
    assert {"hyperscore", "rank"}.issubset(candidates.feature_columns)


def test_exported_proforma_pin_can_be_read_back(tmp_path):
    candidates = read_pin(FIXTURES / "fragpipe_sample.pin")
    output = tmp_path / "roundtrip.pin"

    write_pin(candidates, output)
    roundtrip = read_pin(output)

    columns = ["sequence", "mods", "mod_sites"]
    pd.testing.assert_frame_equal(roundtrip.df[columns], candidates.df[columns])


def test_pin_reader_rejects_labels_other_than_target_or_decoy(tmp_path):
    pin = tmp_path / "invalid-label.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tfilename\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t0\t1\trun\t1.0\t2\tPEPTIDE\tP1\n"
    )

    with pytest.raises(ValueError, match="-1 or 1"):
        read_pin(pin)


def test_pipeline_aligns_observed_charge_features_across_pepxml_runs(tmp_path):
    from optimhc.core.pipeline import Pipeline

    fixture = (FIXTURES / "multi_rank_sample.pep.xml").read_text()
    first = tmp_path / "run_a.pep.xml"
    second = tmp_path / "run_b.pep.xml"
    first.write_text(fixture)
    second.write_text(
        fixture.replace("run_a", "run_b").replace(
            'assumed_charge="2"', 'assumed_charge="3"'
        )
    )
    pipeline = Pipeline(
        {
            "inputType": "pepxml",
            "inputFile": [str(first), str(second)],
            "outputDir": str(tmp_path / "output"),
            "visualization": False,
            "saveModels": False,
        }
    )

    candidates = pipeline.read_input()

    assert {"charge_2", "charge_3"}.issubset(candidates.feature_columns)
    assert candidates.df["charge_2"].tolist() == [1.0, 1.0, 0.0, 0.0]
    assert candidates.df["charge_3"].tolist() == [0.0, 0.0, 1.0, 1.0]
