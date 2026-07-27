from types import SimpleNamespace

import numpy as np

from optimhc.visualization import plot_features


class _Bar:
    def set_hatch(self, value):
        self.hatch = value


class _Axes:
    def barh(self, labels, values, xerr, color):
        self.labels = labels
        self.values = values
        self.errors = xerr
        return [_Bar() for _ in labels]

    def set_xlabel(self, label):
        self.xlabel = label

    def set_ylabel(self, label):
        self.ylabel = label


def test_linear_feature_importance_matches_main_behavior(monkeypatch):
    axes = _Axes()
    monkeypatch.setattr(plot_features.plt, "subplots", lambda **kwargs: (object(), axes))
    monkeypatch.setattr(plot_features, "save_or_show_plot", lambda *args, **kwargs: None)
    models = [
        SimpleNamespace(estimator=SimpleNamespace(coef_=np.array([[-4.0, 2.0], [0.0, -2.0]]))),
        SimpleNamespace(estimator=SimpleNamespace(coef_=np.array([[2.0, 6.0], [2.0, -2.0]]))),
    ]

    plot_features.plot_feature_importance(
        models,
        feature_columns=["feature_a", "feature_b"],
        error=True,
    )

    np.testing.assert_allclose(axes.values, [2.0, 3.0])
    np.testing.assert_allclose(axes.errors, [0.0, 1.0])
