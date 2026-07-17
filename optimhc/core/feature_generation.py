import gc
import logging

import optimhc.feature  # noqa: F401 -- triggers generator registration
from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)

_BINDING_GENERATORS = frozenset({"NetMHCpan", "NetMHCIIpan", "MHCflurry"})


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
    dict
        Raw predictions from configured binding generators, keyed by generator name.
        Empty when ``keepIntermediate`` is disabled.
    """
    feature_generators = config.get("featureGenerator", None)
    if not feature_generators:
        return {}

    keep_intermediate = config.get("keepIntermediate", True)
    raw_predictions = {}

    for generator_config in feature_generators:
        if not isinstance(generator_config, dict):
            logger.warning("Feature generator config is not a dictionary, skipping...")
            continue

        name = generator_config.get("name")
        params = generator_config.get("params", {})

        logger.info(f"Generating features with {name}...")
        generator_cls = feature_generator_factory.get_generator(name)
        generator = generator_cls.from_config(psms, config, params)
        generator.apply(psms, source=name)

        if keep_intermediate and name in _BINDING_GENERATORS:
            predictions = generator.raw_predictions
            if predictions is not None:
                raw_predictions[name] = predictions.copy()

        del generator
        gc.collect()

    return raw_predictions


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
