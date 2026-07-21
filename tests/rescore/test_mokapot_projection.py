import mokapot
import numpy as np
import pandas as pd
import pytest

from optimhc.psm_container import PsmContainer
from optimhc.rescore.model import PercolatorModel
from optimhc.rescore.mokapot import (
    convert_to_mokapot_dataset,
    mokapot_feature_columns,
    rescore,
    to_mokapot_dataframe,
    write_pin,
)


def _candidates():
    return PsmContainer(
        pd.DataFrame(
            {
                "psm_id": [0, 1, 2],
                "run": ["run_a", "run_a", "run_b"],
                "scan": [10, 10, 10],
                "rank": [1, 2, 1],
                "sequence": ["PEPTIDE", "EDITPEP", "PEPTIDE"],
                "mods": ["", "", ""],
                "mod_sites": ["", "", ""],
                "charge": [2, 2, 2],
                "proteins": ["P1", "DECOY_P1", "P1"],
                "is_decoy": [False, True, False],
                "score": [3.0, 1.0, 4.0],
                "retention_time": [12.0, 12.0, 18.0],
                "exp_mass": [900.0, 900.0, 900.0],
            }
        ),
        feature_columns=("score", "rank"),
    )



def test_candidates_for_one_spectrum_share_the_mokapot_spectrum_key():
    projection = to_mokapot_dataframe(_candidates())

    assert projection["SpecId"].tolist() == [
        "run_a.10.10.2",
        "run_a.10.10.2",
        "run_b.10.10.2",
    ]
    assert projection["filename"].tolist() == ["run_a", "run_a", "run_b"]
    assert projection["ScanNr"].tolist() == [10, 10, 10]
    assert projection["Label"].tolist() == [1, -1, 1]
    assert "ret_time" not in projection
    assert "ExpMass" not in projection
    assert projection.columns.tolist() == [
        "SpecId",
        "Label",
        "ScanNr",
        "filename",
        "score",
        "rank",
        "Charge",
        "Peptide",
        "Proteins",
    ]


def test_in_memory_and_serialized_pin_create_equivalent_mokapot_datasets(tmp_path):
    candidates = _candidates()
    in_memory = convert_to_mokapot_dataset(candidates)
    pin_path = tmp_path / "candidates.pin"
    write_pin(candidates, pin_path)
    serialized = mokapot.read_pin(pin_path)

    assert in_memory._spectrum_columns == serialized._spectrum_columns == (
        "filename",
        "ScanNr",
    )
    assert in_memory._feature_columns == serialized._feature_columns
    assert in_memory._feature_columns == ("score", "rank", "Charge")
    pd.testing.assert_frame_equal(in_memory.data, serialized.data)


def test_rescore_passes_explicit_seed_to_mokapot(monkeypatch):
    received = {}

    def fake_brew(dataset, **kwargs):
        received.update(kwargs)
        return "results", "models"

    monkeypatch.setattr("optimhc.rescore.mokapot.mokapot_lib.brew", fake_brew)

    assert rescore(_candidates(), model="model", rng=1) == ("results", "models")
    assert received["model"] == "model"
    assert received["rng"] == 1


def test_rank_cannot_be_removed_from_selected_rescoring_features():
    projection = to_mokapot_dataframe(_candidates(), feature_columns=("score",))

    assert [column for column in ("score", "rank") if column in projection] == [
        "score",
        "rank",
    ]
    assert mokapot_feature_columns(_candidates(), ("score",)) == (
        "score",
        "rank",
        "Charge",
    )


def test_charge_one_hot_features_prevent_scalar_charge_inference():
    candidates = _candidates()
    candidates.add_features(
        pd.DataFrame({"psm_id": [0, 1, 2], "charge_2": [1.0, 1.0, 1.0]}),
        on="psm_id",
        columns=("charge_2",),
    )

    dataset = convert_to_mokapot_dataset(candidates)

    assert dataset._feature_columns == ("score", "rank", "charge_2")
    assert mokapot_feature_columns(candidates) == dataset._feature_columns


def test_charge_one_hot_projection_round_trips_without_scalar_charge(tmp_path):
    candidates = _candidates()
    candidates.add_features(
        pd.DataFrame({"psm_id": [0, 1, 2], "charge_2": [1.0, 1.0, 1.0]}),
        on="psm_id",
        columns=("charge_2",),
    )

    in_memory = convert_to_mokapot_dataset(candidates)
    pin_path = tmp_path / "charge-one-hot.pin"
    projection = write_pin(candidates, pin_path)
    serialized = mokapot.read_pin(pin_path)

    assert "Charge" not in projection
    assert in_memory._feature_columns == serialized._feature_columns
    pd.testing.assert_frame_equal(in_memory.data, serialized.data)


@pytest.mark.parametrize("invalid", ["not-numeric", np.inf, np.nan])
def test_rescoring_boundary_rejects_invalid_declared_features(invalid):
    candidates = _candidates()
    candidates.df["score"] = candidates.df["score"].astype(object)
    candidates.df.loc[0, "score"] = invalid

    with pytest.raises(ValueError, match="finite numeric"):
        to_mokapot_dataframe(candidates)


def test_seeded_rescoring_is_repeatable_and_competes_once_per_spectrum():
    records = []
    for spectrum_index in range(600):
        winning_candidate_is_decoy = spectrum_index % 10 == 0
        for rank in (1, 2):
            is_decoy = (
                winning_candidate_is_decoy if rank == 1 else not winning_candidate_is_decoy
            )
            records.append(
                {
                    "psm_id": len(records),
                    "run": "run_a" if spectrum_index < 300 else "run_b",
                    "scan": spectrum_index % 300 + 1,
                    "rank": rank,
                    "sequence": "PEPTIDEAK" if rank == 1 else "PEPTIDECK",
                    "mods": "",
                    "mod_sites": "",
                    "charge": 2,
                    "proteins": "DECOY_P1" if is_decoy else "P1",
                    "is_decoy": is_decoy,
                    "score": float(20 - rank + spectrum_index / 1000),
                }
            )
    candidates = PsmContainer(
        pd.DataFrame.from_records(records), feature_columns=("score", "rank")
    )
    assert len(convert_to_mokapot_dataset(candidates).data) == 1_200

    confidence = []
    for _ in range(2):
        result, _ = rescore(
            candidates,
            model=PercolatorModel(train_fdr=0.2, n_jobs=1, rng=1),
            test_fdr=0.2,
            rng=1,
        )
        confidence.append(result)

    pd.testing.assert_frame_equal(confidence[0].psms, confidence[1].psms)
    pd.testing.assert_frame_equal(confidence[0].peptides, confidence[1].peptides)
    assert not confidence[0].psms.duplicated(["filename", "ScanNr"]).any()
