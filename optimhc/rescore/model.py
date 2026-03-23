from abc import ABC, abstractmethod

import mokapot.model as mokapot_model
import numpy as np
from mokapot.model import Model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, KFold
from xgboost import XGBClassifier

from optimhc.rescore.factory import rescore_model_factory

# GridSearchCV runtime grows with:
# number of hyperparameter combinations * CV folds * fit cost per model.
# This takes a long time to run
GRID_XGB = {
    "scale_pos_weight": np.logspace(0, 2, 3),
    "max_depth": [3, 5, 7],
    "min_child_weight": [1, 5, 50],
    "gamma": [0, 0.1, 1],
}

GRID_RF = {
    "class_weight": [{0: 1, 1: scale} for scale in np.logspace(0, 2, 3)],
    "max_depth": [3, 5, 7],
    "min_samples_split": [2, 5, 50],
    "min_impurity_decrease": [0, 0.1, 1],
}


class MokapotModelWrapper(Model, ABC):
    """Base class for mokapot-compatible models created from pipeline config."""

    @classmethod
    @abstractmethod
    def from_config(cls, config: dict) -> "MokapotModelWrapper":
        """Build a model instance from full pipeline config."""
        raise NotImplementedError

    @staticmethod
    def _rescore_config(config: dict) -> dict:
        return config.get("rescore", {})


class PercolatorModel(mokapot_model.PercolatorModel, MokapotModelWrapper):
    """Wrapper around mokapot's built-in Percolator model."""

    @classmethod
    def from_config(cls, config: dict) -> "PercolatorModel":
        rescore_cfg = cls._rescore_config(config)
        return cls(
            train_fdr=rescore_cfg.get("trainFDR", 0.01),
            n_jobs=rescore_cfg.get("numJobs", 1),
        )


class XGBoostModel(MokapotModelWrapper):
    def __init__(
        self,
        scaler=None,
        train_fdr=0.01,
        max_iter=10,
        direction=None,
        override=False,
        n_jobs=1,
        rng=None,
    ):
        self.n_jobs = n_jobs
        rng_instance = np.random.default_rng(rng)
        estimator = GridSearchCV(
            # keep estimator single-threaded; GridSearchCV handles parallelism.
            XGBClassifier(random_state=42, n_jobs=1),
            param_grid=GRID_XGB,
            refit=False,
            cv=KFold(3, shuffle=True, random_state=rng_instance.integers(1, 1e6)),
            n_jobs=n_jobs,
            scoring="roc_auc",
        )
        super().__init__(
            estimator=estimator,
            scaler=scaler,
            train_fdr=train_fdr,
            max_iter=max_iter,
            direction=direction,
            override=override,
            rng=rng,
        )

    @classmethod
    def from_config(cls, config: dict) -> "XGBoostModel":
        rescore_cfg = cls._rescore_config(config)
        return cls(
            train_fdr=rescore_cfg.get("trainFDR", 0.01),
            n_jobs=rescore_cfg.get("numJobs", 1),
        )


class RandomForestModel(MokapotModelWrapper):
    def __init__(
        self,
        scaler=None,
        train_fdr=0.01,
        max_iter=10,
        direction=None,
        override=False,
        n_jobs=1,
        rng=None,
    ):
        self.n_jobs = n_jobs
        rng_instance = np.random.default_rng(rng)
        estimator = GridSearchCV(
            # keep estimator single-threaded; GridSearchCV handles parallelism.
            RandomForestClassifier(random_state=42, n_jobs=1),
            param_grid=GRID_RF,
            refit=False,
            cv=KFold(3, shuffle=True, random_state=rng_instance.integers(1, 1e6)),
            n_jobs=n_jobs,
            scoring="roc_auc",
        )
        super().__init__(
            estimator=estimator,
            scaler=scaler,
            train_fdr=train_fdr,
            max_iter=max_iter,
            direction=direction,
            override=override,
            rng=rng,
        )

    @classmethod
    def from_config(cls, config: dict) -> "RandomForestModel":
        rescore_cfg = cls._rescore_config(config)
        return cls(
            train_fdr=rescore_cfg.get("trainFDR", 0.01),
            n_jobs=rescore_cfg.get("numJobs", 1),
        )


rescore_model_factory.register_model("Percolator", PercolatorModel)
rescore_model_factory.register_model("XGBoost", XGBoostModel)
rescore_model_factory.register_model("RandomForest", RandomForestModel)
