from types import SimpleNamespace

import pandas as pd
import pytest

from optimhc.core.pipeline import Pipeline
from optimhc.output.flashlfq import format_flashlfq, write_flashlfq
from optimhc.psm_container import PsmContainer


def _candidates(include_retention_time=True, include_calc_mass=True):
    data = {
        "psm_id": [0, 1],
        "run": ["run_a", "run_b"],
        "scan": [10, 20],
        "rank": [1, 1],
        "sequence": ["PEPTM", "OTHER"],
        "mods": ["Oxidation@M", ""],
        "mod_sites": ["5", ""],
        "charge": [2, 3],
        "proteins": ["P1;P2", "P3"],
        "is_decoy": [False, False],
        "score": [3.0, 2.0],
    }
    if include_retention_time:
        data["retention_time"] = [120.0, 180.0]
    if include_calc_mass:
        data["calc_mass"] = [600.25, 700.5]
    return PsmContainer(pd.DataFrame(data), feature_columns=("score", "rank"))


def _peptide_results():
    return pd.DataFrame(
        {
            "SpecId": ["run_a.10.10.2", "run_b.20.20.3"],
            "Peptide": ["PEPTM[UNIMOD:35]", "OTHER"],
            "mokapot q-value": [0.005, 0.02],
        }
    )


def test_flashlfq_formats_accepted_peptides_from_psm_dataframe():
    output = format_flashlfq(_candidates(), _peptide_results(), fdr=0.01)

    assert output.to_dict("records") == [
        {
            "File Name": "run_a",
            "Base Sequence": "PEPTM",
            "Full Sequence": "PEPTM[UNIMOD:35]",
            "Peptide Monoisotopic Mass": 600.25,
            "Scan Retention Time": 2.0,
            "Precursor Charge": 2,
            "Protein Accession": "P1;P2",
        }
    ]


def test_flashlfq_can_derive_missing_calculated_mass():
    output = format_flashlfq(_candidates(include_calc_mass=False), _peptide_results(), fdr=0.01)

    assert output["Peptide Monoisotopic Mass"].iat[0] == pytest.approx(589.2418, abs=1e-3)


def test_flashlfq_requires_retention_time():
    with pytest.raises(ValueError, match="retention_time"):
        format_flashlfq(_candidates(include_retention_time=False), _peptide_results(), fdr=0.01)


def test_flashlfq_writer_accepts_mokapot_confidence_object(tmp_path):
    results = SimpleNamespace(peptides=_peptide_results())
    path = tmp_path / "optimhc.FlashLFQ.txt"

    written = write_flashlfq(_candidates(), results, path, fdr=0.01)

    assert written == path
    assert pd.read_csv(path, sep="\t").shape == (1, 7)


@pytest.mark.parametrize("enabled, expected_calls", [(True, 1), (False, 0)])
def test_pipeline_honors_to_flashlfq_config(tmp_path, monkeypatch, enabled, expected_calls):
    input_path = tmp_path / "input.pin"
    input_path.touch()
    pipeline = Pipeline(
        {
            "inputType": "pin",
            "inputFile": [str(input_path)],
            "outputDir": str(tmp_path / "output"),
            "saveModels": False,
            "toFlashLFQ": enabled,
        }
    )
    calls = []

    monkeypatch.setattr("optimhc.rescore.mokapot.write_pin", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        "optimhc.output.flashlfq.write_flashlfq",
        lambda *args, **kwargs: calls.append((args, kwargs)),
    )
    results = SimpleNamespace(to_txt=lambda **kwargs: None)

    pipeline.save_results(object(), results, [])

    assert len(calls) == expected_calls
    if calls:
        assert calls[0][1]["fdr"] == 0.01
        assert str(calls[0][0][2]).endswith("optimhc.FlashLFQ.txt")


def test_pipeline_rejects_missing_retention_time_before_feature_generation(tmp_path, monkeypatch):
    input_path = tmp_path / "input.pin"
    input_path.touch()
    pipeline = Pipeline(
        {
            "inputType": "pin",
            "inputFile": [str(input_path)],
            "outputDir": str(tmp_path / "output"),
            "saveModels": False,
        }
    )
    monkeypatch.setattr(
        pipeline,
        "read_input",
        lambda: _candidates(include_retention_time=False),
    )
    monkeypatch.setattr(
        pipeline,
        "_generate_features",
        lambda psms: pytest.fail("Feature generation should not run without retention time."),
    )

    with pytest.raises(ValueError, match="FlashLFQ output requires retention time"):
        pipeline.run()
