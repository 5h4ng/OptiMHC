"""Visualizations for explicitly declared rescoring features."""

import logging

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from optimhc.psm_container import PsmContainer
from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)


def plot_feature_importance(
    models,
    feature_columns,
    save_path=None,
    sort=False,
    error=False,
    **kwargs,
):
    """Plot average feature importance across trained folds."""
    features = list(feature_columns)
    estimator = models[0].estimator
    if hasattr(estimator, "coef_"):
        values = np.asarray([model.estimator.coef_.mean(axis=0) for model in models])
        importance = np.abs(values).mean(axis=0)
        signs = values.mean(axis=0)
    elif hasattr(estimator, "feature_importances_"):
        values = np.asarray([model.estimator.feature_importances_ for model in models])
        importance = values.mean(axis=0)
        signs = np.ones_like(importance)
    else:
        raise ValueError("Model does not expose coefficients or feature importances.")
    if len(features) != len(importance):
        raise ValueError("Feature names do not match the trained model dimensions.")

    order = np.argsort(-importance) if sort else np.arange(len(features))
    errors = values.std(axis=0)[order] if error else None
    fig, ax = plt.subplots(
        figsize=kwargs.get("figsize", (15, 10)),
        dpi=kwargs.get("dpi", 300),
    )
    bars = ax.barh(
        [features[index] for index in order],
        importance[order],
        xerr=errors,
        color=sns.color_palette(kwargs.get("palette", "Set2"), n_colors=1)[0],
    )
    for bar, sign in zip(bars, signs[order]):
        if sign < 0:
            bar.set_hatch("\\\\")
    ax.set_xlabel("Average Feature Importance")
    ax.set_ylabel("Feature")
    save_or_show_plot(save_path, logger)


def visualize_feature_correlation(psms: PsmContainer, save_path=None, **kwargs):
    """Visualize correlations among the declared model features."""
    features = list(psms.feature_columns)
    count = len(features)
    correlation = psms.df[features].corr().stack().reset_index(name="correlation")
    plot = sns.relplot(
        data=correlation,
        x="level_0",
        y="level_1",
        hue="correlation",
        size="correlation",
        palette="vlag",
        hue_norm=(-1, 1),
        edgecolor=".7",
        height=kwargs.get("height", max(8, min(20, 8 + count * 0.15))),
        sizes=(50, 250),
        size_norm=(-0.2, 0.8),
    )
    plot.set(xlabel="", ylabel="", aspect="equal")
    plot.despine(left=True, bottom=True)
    for label in plot.ax.get_xticklabels():
        label.set_rotation(90)
    plot.figure.dpi = kwargs.get("dpi", 300)
    plt.tight_layout()
    save_or_show_plot(save_path, logger)
