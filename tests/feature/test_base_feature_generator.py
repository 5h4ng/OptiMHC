import pandas as pd
import pytest

from optimhc.core.feature_generation import generate_features, select_feature_groups
from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.psm_container import PsmContainer


class SequenceFeature(BaseFeatureGenerator):
    @property
    def feature_columns(self):
        return ["score"]

    @property
    def id_column(self):
        return ["Peptide"]

    def generate_features(self):
        return pd.DataFrame(
            {"Peptide": ["PEPTIDE", "PEPTIDE"], "score": [1.0, 1.0]}
        )


def test_default_sequence_feature_join_collapses_identical_duplicate_predictions():
    candidates = PsmContainer(
        pd.DataFrame(
            {
                "psm_id": [0, 1],
                "run": ["run", "run"],
                "scan": [1, 2],
                "rank": [1, 1],
                "sequence": ["PEPTIDE", "PEPTIDE"],
                "mods": ["", ""],
                "mod_sites": ["", ""],
                "charge": [2, 2],
                "proteins": ["P1", "P1"],
                "is_decoy": [False, False],
            }
        )
    )

    SequenceFeature().apply(candidates)

    assert candidates.df["score"].tolist() == [1.0, 1.0]


def test_generator_declares_its_feature_group():
    assert SequenceFeature().feature_groups("Basic") == {"Basic": ("score",)}


def test_feature_groups_resolve_to_columns_in_source_order():
    manifest = {
        "Original": ("xcorr", "rank"),
        "Basic": ("length", "entropy"),
    }

    assert select_feature_groups(manifest, ["Basic", "Original"]) == (
        "length",
        "entropy",
        "xcorr",
        "rank",
    )


def test_unknown_feature_group_fails_explicitly():
    manifest = {"Original": ("xcorr", "rank")}

    with pytest.raises(ValueError, match="Unknown feature source.*Basic"):
        select_feature_groups(manifest, ["Basic"])


def test_feature_generation_returns_run_local_manifest():
    candidates = PsmContainer(
        pd.DataFrame(
            {
                "psm_id": [0, 1],
                "run": ["run", "run"],
                "scan": [1, 2],
                "rank": [1, 1],
                "sequence": ["PEPTIDE", "SEQUENCE"],
                "mods": ["", ""],
                "mod_sites": ["", ""],
                "charge": [2, 2],
                "proteins": ["P1", "P2"],
                "is_decoy": [False, True],
            }
        )
    )

    generated = generate_features(
        candidates,
        {
            "featureGenerator": [{"name": "Basic"}],
            "removePreNxtAA": False,
            "keepIntermediate": False,
        },
    )

    assert generated.feature_groups["Original"] == ("rank",)
    assert generated.feature_groups["Basic"] == tuple(
        column for column in candidates.feature_columns if column != "rank"
    )
