import gc
import logging

import optimhc.feature  # noqa: F401 -- triggers generator registration
from optimhc.feature.factory import feature_generator_factory

logger = logging.getLogger(__name__)


def generate_features(psms, config):
    """
    Generate features from different generators according to the configuration.

    Parameters
    ----------
    psms : PsmContainer
        A container object holding PSMs and relevant data.
    config : dict
        Configuration dictionary loaded from YAML or CLI.
    """
    feature_generators = config.get("featureGenerator", None)
    if not feature_generators:
        return

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
        gc.collect()
