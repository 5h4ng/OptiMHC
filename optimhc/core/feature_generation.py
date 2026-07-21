import gc
import importlib
import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import pandas as pd

from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)

_BINDING_GENERATORS = frozenset({"NetMHCpan", "NetMHCIIpan", "MHCflurry"})
_GENERATOR_MODULES = {
    "Basic": "basic",
    "DeepLC": "deeplc",
    "MHCflurry": "mhcflurry",
    "NetMHCIIpan": "netmhciipan",
    "NetMHCpan": "netmhcpan",
    "OverlappingPeptide": "overlapping_peptide",
    "PWM": "pwm",
    "SpectralSimilarity": "spectral_similarity",
}


@dataclass(frozen=True)
class FeatureGenerationResult:
    """Run-local feature groups plus optional raw predictor outputs."""

    feature_groups: dict[str, tuple[str, ...]]
    raw_predictions: dict[str, pd.DataFrame]


def select_feature_groups(
    feature_groups: Mapping[str, Sequence[str]],
    sources: Sequence[str],
) -> tuple[str, ...]:
    """Resolve legacy experiment source names to explicit feature columns."""
    unknown = [source for source in sources if source not in feature_groups]
    if unknown:
        available = ", ".join(feature_groups)
        raise ValueError(f"Unknown feature source(s): {unknown}. Available sources: {available}")

    selected = []
    seen = set()
    for source in sources:
        for column in feature_groups[source]:
            if column not in seen:
                selected.append(column)
                seen.add(column)
    return tuple(selected)


def generate_features(psms, config):
    """
    Generate features from different generators according to the configuration.

    Parameters
    ----------
    psms : PsmContainer
        A container object holding PSMs and relevant data.
    config : dict
        Configuration dictionary loaded from YAML or CLI.

    Returns
    -------
    FeatureGenerationResult
        Generator-to-column groups for experiment selection and optional raw
        binding predictions for intermediate output.
    """
    feature_generators = config.get("featureGenerator", None)
    feature_groups = {"Original": tuple(psms.feature_columns)}
    raw_predictions = {}
    if not feature_generators:
        return FeatureGenerationResult(feature_groups, raw_predictions)

    keep_intermediate = config.get("keepIntermediate", True)

    for generator_config in feature_generators:
        if not isinstance(generator_config, dict):
            logger.warning("Feature generator config is not a dictionary, skipping...")
            continue

        name = generator_config.get("name")
        params = generator_config.get("params", {})

        logger.info(f"Generating features with {name}...")
        try:
            module = _GENERATOR_MODULES[name]
        except KeyError as error:
            raise ValueError(f"Unknown feature generator: '{name}'.") from error
        importlib.import_module(f"optimhc.feature.{module}")
        generator_cls = feature_generator_factory.get_generator(name)
        generator = generator_cls.from_config(psms, config, params)
        previous_columns = set(psms.feature_columns)
        generator.apply(psms)

        generated_columns = tuple(
            column for column in psms.feature_columns if column not in previous_columns
        )
        declared_groups = generator.feature_groups(name)
        declared_columns = tuple(
            column for columns in declared_groups.values() for column in columns
        )
        if set(generated_columns) != set(declared_columns):
            raise ValueError(
                f"Generator '{name}' added {generated_columns}, but declared {declared_columns}."
            )
        duplicate_groups = set(feature_groups).intersection(declared_groups)
        if duplicate_groups:
            raise ValueError(f"Feature groups declared more than once: {duplicate_groups}")
        feature_groups.update(
            {group: tuple(columns) for group, columns in declared_groups.items()}
        )

        if keep_intermediate and name in _BINDING_GENERATORS:
            predictions = generator.raw_predictions
            if predictions is not None:
                raw_predictions[name] = predictions.copy()

        del generator
        gc.collect()

    return FeatureGenerationResult(feature_groups, raw_predictions)


def _build_ba_parquet(raw_predictions, output_path):
    """Assemble BA.parquet from collected raw binding predictions.

    Parameters
    ----------
    raw_predictions : dict
        Dict mapping generator name to raw prediction DataFrame.
    output_path : str
        Path to write the parquet file.
    """
    frames = {}

    for name, df in raw_predictions.items():
        df = df.copy()
        if name in ("NetMHCpan", "NetMHCIIpan"):
            prefix = name.lower()
            valid = df.dropna(subset=["percentile_rank"])
            if valid.empty:
                logger.warning("No valid predictions in %s, skipping.", name)
                continue
            idx = valid.groupby("peptide")["percentile_rank"].idxmin()
            best = valid.loc[idx].copy()
            best = best.rename(
                columns={
                    "allele": f"{prefix}_allele",
                    "affinity": f"{prefix}_affinity",
                    "percentile_rank": f"{prefix}_ba_rank",
                }
            )
            frames[name] = best[
                ["peptide", f"{prefix}_allele", f"{prefix}_affinity", f"{prefix}_ba_rank"]
            ]

        elif name == "MHCflurry":
            mhc = df.rename(
                columns={
                    "best_allele": "mhcflurry_allele",
                    "affinity": "mhcflurry_affinity",
                    "presentation_percentile": "mhcflurry_el_rank",
                }
            )
            cols = ["peptide", "mhcflurry_allele", "mhcflurry_affinity", "mhcflurry_el_rank"]
            frames[name] = mhc[[c for c in cols if c in mhc.columns]]

    if not frames:
        logger.warning("No binding predictions to save.")
        return

    frame_list = list(frames.values())
    ba = frame_list[0]
    for other in frame_list[1:]:
        ba = ba.merge(other, on="peptide", how="outer")

    ba.set_index("peptide", inplace=True)
    ba.to_parquet(output_path)
    logger.info(
        "BA.parquet saved to %s (%d peptides, %d columns).",
        output_path,
        len(ba),
        len(ba.columns),
    )
