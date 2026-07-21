"""Feature generation interfaces.

Concrete predictors are imported only when selected in pipeline configuration.
"""

from optimhc.feature.base_feature_generator import BaseFeatureGenerator
from optimhc.feature.factory import feature_generator_factory

__all__ = ["BaseFeatureGenerator", "feature_generator_factory"]
