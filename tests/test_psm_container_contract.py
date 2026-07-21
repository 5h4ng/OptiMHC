import pandas as pd
import pytest

from optimhc.psm_container import PsmContainer


def test_container_exposes_the_psm_dataframe_and_features():
    candidates = pd.DataFrame(
        {
            "psm_id": [0, 1],
            "run": ["run_a", "run_a"],
            "scan": [10, 10],
            "rank": [1, 2],
            "sequence": ["PEPTIDE", "EDITPEP"],
            "mods": ["", ""],
            "mod_sites": ["", ""],
            "charge": [2, 2],
            "proteins": ["P1", "DECOY_P1"],
            "is_decoy": [False, True],
            "search_score": [3.0, 1.0],
        }
    )

    container = PsmContainer(candidates, feature_columns=("search_score", "rank"))

    assert container.df is candidates
    assert container.feature_columns == ("search_score", "rank")
    assert len(container) == 2


def test_feature_author_can_attach_sequence_evidence_explicitly():
    candidates = pd.DataFrame(
        {
            "psm_id": [0, 1, 2],
            "run": ["run_a", "run_a", "run_b"],
            "scan": [10, 11, 10],
            "rank": [1, 1, 1],
            "sequence": ["PEPTIDE", "EDITPEP", "PEPTIDE"],
            "mods": ["", "", ""],
            "mod_sites": ["", "", ""],
            "charge": [2, 2, 3],
            "proteins": ["P1", "P2", "P1"],
            "is_decoy": [False, True, False],
        }
    )
    container = PsmContainer(candidates)

    container.add_features(
        pd.DataFrame(
            {
                "sequence": ["PEPTIDE", "EDITPEP"],
                "binding_score": [0.1, 0.9],
            }
        ),
        on="sequence",
        columns=("binding_score",),
    )

    assert container.df["binding_score"].tolist() == [0.1, 0.9, 0.1]
    assert container.feature_columns == ("rank", "binding_score")


def test_feature_join_rejects_incomplete_key_coverage():
    candidates = pd.DataFrame(
        {
            "psm_id": [0, 1],
            "run": ["run_a", "run_a"],
            "scan": [10, 11],
            "rank": [1, 1],
            "sequence": ["PEPTIDE", "EDITPEP"],
            "mods": ["", ""],
            "mod_sites": ["", ""],
            "charge": [2, 2],
            "proteins": ["P1", "P2"],
            "is_decoy": [False, True],
        }
    )
    container = PsmContainer(candidates)

    with pytest.raises(ValueError, match="exactly cover"):
        container.add_features(
            pd.DataFrame({"sequence": ["PEPTIDE"], "score": [0.1]}),
            on="sequence",
            columns=("score",),
        )

    assert "score" not in container.df


def test_feature_join_stores_normalized_numeric_values():
    candidates = pd.DataFrame(
        {
            "psm_id": [0],
            "run": ["run_a"],
            "scan": [10],
            "rank": [1],
            "sequence": ["PEPTIDE"],
            "mods": [""],
            "mod_sites": [""],
            "charge": [2],
            "proteins": ["P1"],
            "is_decoy": [False],
        }
    )
    container = PsmContainer(candidates)

    container.add_features(
        pd.DataFrame({"psm_id": [0], "score": ["1.5"]}),
        on="psm_id",
        columns=("score",),
    )

    assert container.df["score"].tolist() == [1.5]
    assert pd.api.types.is_numeric_dtype(container.df["score"])
