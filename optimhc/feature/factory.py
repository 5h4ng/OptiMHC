import logging
from typing import Dict, List, Type

from optimhc.feature.base_feature_generator import BaseFeatureGenerator

logger = logging.getLogger(__name__)


class FeatureGeneratorFactory:
    """Registry and factory for feature generators.

    Each generator module registers itself at import time by calling
    ``feature_generator_factory.register_generator(name, cls)``.
    The orchestrator retrieves generators by name via ``get_generator(name)``.
    """

    def __init__(self):
        self._registry: Dict[str, Type[BaseFeatureGenerator]] = {}

    def register_generator(self, name: str, generator_class: Type[BaseFeatureGenerator]):
        """Register a feature generator class under *name*."""
        self._registry[name] = generator_class

    def get_generator(self, name: str) -> Type[BaseFeatureGenerator]:
        """Return the generator class registered under *name*."""
        if name not in self._registry:
            raise ValueError(
                f"Unknown feature generator: '{name}'. Available: {sorted(self._registry.keys())}"
            )
        return self._registry[name]

    def list_generators(self) -> List[str]:
        """Return sorted list of registered generator names."""
        return sorted(self._registry.keys())


feature_generator_factory = FeatureGeneratorFactory()
