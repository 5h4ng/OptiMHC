from abc import ABC, abstractmethod
from typing import List

import pandas as pd

from optimhc.psm_container import PsmContainer


class BaseFeatureGenerator(ABC):
    """Abstract base class for all feature generators in the rescoring pipeline.

    Subclasses must implement:
    - ``feature_columns`` -- names of generated feature columns
    - ``id_column`` -- merge key column(s)
    - ``generate_features()`` -- pure computation, returns a DataFrame
    - ``from_config()`` -- construct an instance from pipeline config

    The default ``apply()`` merges features by peptide column.
    Override it for index-based merges, composite keys, or post-processing.
    """

    @property
    @abstractmethod
    def feature_columns(self) -> List[str]:
        """Return a list of feature column names produced by this generator."""
        ...

    @property
    @abstractmethod
    def id_column(self) -> List[str]:
        """Return the column(s) used as merge key(s) with the PsmContainer."""
        ...

    @abstractmethod
    def generate_features(self) -> pd.DataFrame:
        """Generate features and return them as a DataFrame."""
        ...

    @classmethod
    def from_config(
        cls,
        psms: PsmContainer,
        config: dict,
        params: dict,
    ) -> "BaseFeatureGenerator":
        """Construct a generator instance from pipeline configuration.

        Parameters
        ----------
        psms : PsmContainer
            The PSM container with all current data.
        config : dict
            The full pipeline configuration.
        params : dict
            Generator-specific parameters from
            ``config["featureGenerator"][i]["params"]``.
        """
        raise NotImplementedError(f"{cls.__name__} must implement from_config()")

    def apply(self, psms: PsmContainer, source: str) -> None:
        """Generate features and merge them into the PsmContainer.

        The default implementation merges by peptide column using
        ``add_features()``.  Override for different merge strategies
        (index-based, composite key) or additional post-processing.

        Parameters
        ----------
        psms : PsmContainer
            The PSM container to add features to (modified in-place).
        source : str
            Name of this feature source (e.g. ``"Basic"``, ``"PWM"``).
        """
        features = self.generate_features()
        psms.add_features(
            features,
            psms_key=psms.peptide_column,
            feature_key=self.id_column,
            source=source,
        )
