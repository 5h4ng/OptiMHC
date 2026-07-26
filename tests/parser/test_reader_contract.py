from pathlib import Path

import pandas as pd
import pytest

from optimhc.parser import read_pepxml, read_pin
from optimhc.parser.modifications import _select_best_match, modification_from_delta
from optimhc.rescore.mokapot import write_pin

FIXTURES = Path(__file__).parent


def test_pin_reader_returns_normalized_candidates():
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
    assert candidates.df["matched_ions_ratio"].tolist() == [0.5, 0.25]
    assert {"hyperscore", "rank"}.issubset(candidates.feature_columns)


def test_pepxml_reader_skips_matched_ions_ratio_for_nonpositive_totals(tmp_path, caplog):
    fixture = (FIXTURES / "multi_rank_sample.pep.xml").read_text()
    path = tmp_path / "zero-total-ions.pep.xml"
    path.write_text(fixture.replace('tot_num_ions="8"', 'tot_num_ions="0"', 1))

    candidates = read_pepxml(path)

    assert "matched_ions_ratio" not in candidates.df
    assert "matched_ions_ratio" not in candidates.feature_columns
    assert "Skipping matched_ions_ratio" in caplog.text


def test_exported_proforma_pin_can_be_read_back(tmp_path):
    candidates = read_pin(FIXTURES / "fragpipe_sample.pin")
    output = tmp_path / "roundtrip.pin"

    write_pin(candidates, output)
    roundtrip = read_pin(output)

    columns = ["sequence", "mods", "mod_sites"]
    pd.testing.assert_frame_equal(roundtrip.df[columns], candidates.df[columns])


@pytest.mark.parametrize(
    ("peptide", "sequence", "mods", "mod_sites"),
    [
        ("K.PEPTIDE.R", "PEPTIDE", "", ""),
        ("PEPTIDE", "PEPTIDE", "", ""),
        ("K.PEPM[15.9949]K.R", "PEPMK", "Oxidation@M", "4"),
        ("PEPM[+15.9949]K", "PEPMK", "Oxidation@M", "4"),
        (
            "K.Q[-17.0265]PEPTIDE.R",
            "QPEPTIDE",
            "Gln->pyro-Glu@Q^Any_N-term",
            "0",
        ),
        (
            "K.C[-17.0265]PEPTIDE.R",
            "CPEPTIDE",
            "Ammonia-loss@C^Any_N-term",
            "0",
        ),
        (
            "F.VTVQGRAIC[119.0041]SDPNNKRVKN4.A",
            "VTVQGRAICSDPNNKRVKN",
            "Cysteinyl@C",
            "9",
        ),
        (
            "K.n[42.0106]PEPM[15.9949]TIDEc[-0.9840].R",
            "PEPMTIDE",
            "Acetyl@Any_N-term;Oxidation@M;Amidated@Any_C-term",
            "0;4;-1",
        ),
        (
            "[+42.0106]-PEPM[+15.9949]TIDE-[-0.9840]",
            "PEPMTIDE",
            "Acetyl@Any_N-term;Oxidation@M;Amidated@Any_C-term",
            "0;4;-1",
        ),
    ],
)
def test_pin_reader_supports_numeric_peptide_dialects(
    tmp_path, peptide, sequence, mods, mod_sites
):
    pin = tmp_path / "peptide-dialect.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tscore\tCharge\tPeptide\tProteins\n"
        f"run.1.1.2\t1\t1\t5.0\t2\t{peptide}\tP1\n"
    )

    candidate = read_pin(pin).df.iloc[0]

    assert candidate["sequence"] == sequence
    assert candidate["mods"] == mods
    assert candidate["mod_sites"] == mod_sites


def test_pin_reader_rejects_symbol_modifications(tmp_path):
    pin = tmp_path / "symbol-modification.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t1\t1\t5.0\t2\tK.PEPM*K.R\tP1\n"
    )

    with pytest.raises(ValueError, match="Unsupported PIN peptide syntax"):
        read_pin(pin)


def test_pin_reader_rejects_terminal_mass_on_wrong_residue(tmp_path):
    pin = tmp_path / "wrong-terminal-residue.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t1\t1\t5.0\t2\tK.M[-17.0265]PEPTIDE.R\tP1\n"
    )

    with pytest.raises(ValueError, match="Unknown modification mass"):
        read_pin(pin)


def test_ambiguous_modification_mapping_logs_once(caplog):
    modification_from_delta.cache_clear()

    first = modification_from_delta(-15.994915, residue="S")
    second = modification_from_delta(-15.994915, residue="S")

    assert first == second == "Deoxy@S"
    assert caplog.text.count("Multiple modification definitions match") == 1
    assert "Deoxy@S, Ser->Ala@S" in caplog.text


def test_terminal_mapping_prefers_any_over_protein_terminus():
    modification_from_delta.cache_clear()

    result = modification_from_delta(42.010565, residue="Q", site=0)

    assert result == "Acetyl@Any_N-term"


def test_modification_mapping_prefers_smallest_mass_error():
    candidates = pd.DataFrame(
        {
            "mod_name": ["First@M", "Closer@M"],
            "unimod_mass": [15.0, 15.003],
        }
    )

    result = _select_best_match(candidates, delta_mass=15.0029)

    assert result == "Closer@M"


def test_pin_reader_rejects_labels_other_than_target_or_decoy(tmp_path):
    pin = tmp_path / "invalid-label.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tfilename\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t0\t1\trun\t1.0\t2\tPEPTIDE\tP1\n"
    )

    with pytest.raises(ValueError, match="-1 or 1"):
        read_pin(pin)


def test_pin_reader_ignores_unknown_nonnumeric_columns(tmp_path, caplog):
    pin = tmp_path / "extra-metadata.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tSearchEngine\tscore\tCharge\tPeptide\tProteins\n"
        "run.1.1.2\t1\t1\tComet\t5.0\t2\tK.PEPTIDE.R\tP1\n"
    )

    candidates = read_pin(pin)

    assert "SearchEngine" not in candidates.df
    assert "SearchEngine" not in candidates.feature_columns
    assert "score" in candidates.feature_columns
    assert "Ignoring non-numeric PIN column 'SearchEngine'" in caplog.text


@pytest.mark.parametrize(
    ("filename", "expected_run"),
    [
        ("run_a.mzML", "run_a"),
        ("run_a.raw", "run_a.raw"),
        ("run_a.RAW", "run_a.RAW"),
    ],
)
def test_pin_reader_recognizes_rt_and_filename_metadata(
    tmp_path,
    filename,
    expected_run,
):
    pin = tmp_path / "percolator-metadata.pin"
    pin.write_text(
        "SpecId\tLabel\tScanNr\tfilename\trt\tscore\tCharge\tPeptide\tProteins\n"
        f"psm-1\t1\t1\t{filename}\t12.0\t5.0\t2\tK.PEPTIDE.R\tP1\n"
    )

    candidates = read_pin(pin)

    assert candidates.df.loc[0, "run"] == expected_run
    assert candidates.df.loc[0, "retention_time"] == 720.0
    assert "filename" not in candidates.feature_columns
    assert "rt" not in candidates.feature_columns
    assert "Charge" in candidates.feature_columns


@pytest.mark.parametrize(
    ("retention_times", "expected"),
    [
        ([12.0, 20.0], [720.0, 1200.0]),
        ([12.0, 500.0], [12.0, 500.0]),
    ],
)
def test_pin_reader_uses_500_as_retention_time_threshold(
    tmp_path,
    retention_times,
    expected,
):
    pin = tmp_path / "retention-time.pin"
    rows = [
        "SpecId\tLabel\tScanNr\trt\tscore\tCharge\tPeptide\tProteins",
        (f"run.1.1.2\t1\t1\t{retention_times[0]}\t5.0\t2\tK.PEPTIDE.R\tP1"),
        (f"run.2.2.2\t1\t2\t{retention_times[1]}\t4.0\t2\tK.PEPTIDE.R\tP1"),
    ]
    pin.write_text("\n".join(rows) + "\n")

    candidates = read_pin(pin)

    assert candidates.df["retention_time"].tolist() == expected


def test_pipeline_aligns_observed_charge_features_across_pepxml_runs(tmp_path):
    from optimhc.core.pipeline import Pipeline

    fixture = (FIXTURES / "multi_rank_sample.pep.xml").read_text()
    first = tmp_path / "run_a.pep.xml"
    second = tmp_path / "run_b.pep.xml"
    first.write_text(fixture)
    second.write_text(
        fixture.replace("run_a", "run_b").replace('assumed_charge="2"', 'assumed_charge="3"')
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
