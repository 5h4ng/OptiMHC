from pathlib import Path

import pytest

from optimhc.parser import read_pepxml, read_pin

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


def test_pin_reader_rejects_labels_other_than_target_or_decoy(tmp_path):
    pin = tmp_path / "invalid-label.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tfilename\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t0\t1\trun\t1.0\t2\tPEPTIDE\tP1\n"
    )

    with pytest.raises(ValueError, match="-1 or 1"):
        read_pin(pin)
