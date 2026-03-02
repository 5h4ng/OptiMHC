# rescore/model.py

import numpy as np
from mokapot.model import Model
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, KFold
from xgboost import XGBClassifier

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


class XGBoostPercolatorModel(Model):
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
        grid = GRID_XGB
        estimator = GridSearchCV(
            XGBClassifier(random_state=42),
            param_grid=grid,
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


# TODO: not tested
class RandomForestPercolatorModel(Model):
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
            RandomForestClassifier(random_state=42),
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
