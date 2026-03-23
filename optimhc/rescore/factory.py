from typing import Dict, List, Type

from mokapot.model import Model


class RescoreModelFactory:
    """Registry and factory for rescore models."""

    def __init__(self):
        self._registry: Dict[str, Type[Model]] = {}

    def register_model(self, name: str, model_class: Type[Model]) -> None:
        """Register a model class under *name*."""
        self._registry[name] = model_class

    def get_model(self, name: str) -> Type[Model]:
        """Return the model class registered under *name*."""
        if name not in self._registry:
            raise ValueError(
                f"Unknown rescore model: '{name}'. Available: {sorted(self._registry.keys())}"
            )
        return self._registry[name]

    def list_models(self) -> List[str]:
        """Return sorted list of registered model names."""
        return sorted(self._registry.keys())


rescore_model_factory = RescoreModelFactory()
