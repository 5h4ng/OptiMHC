from typing import Any, Dict, List, Type


class RescoreModelFactory:
    """Registry and factory for rescore models."""

    def __init__(self):
        self._registry: Dict[str, Type[Any]] = {}

    def register_model(self, name: str, model_class: Type[Any]) -> None:
        """Register a model class under *name*."""
        self._registry[name] = model_class

    def get_model(self, name: str) -> Type[Any]:
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
