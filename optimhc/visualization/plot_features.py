import logging

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Patch

from optimhc.psm_container import PsmContainer
from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)


def plot_feature_importance(
    models, rescoring_features, save_path=None, sort=False, error=False, **kwargs
):
    """
    Unified function to plot average feature importance across multiple models.

    This function supports:
      - Linear models (e.g., Linear SVR) which provide an 'estimator' attribute with a 'coef_'.
        The absolute value of the coefficients is used for importance, and hatch patterns are applied
        to differentiate between positive and negative coefficients.
      - XGBoost models which provide a 'feature_importances_' attribute. Since these values are
        always positive, no hatch patterns are applied.

    Parameters
    ----------
    models : list
        A list of model objects.
        For linear models, each model should have an 'estimator' with 'coef_'.
        For XGBoost models, each model should have a 'feature_importances_' attribute.
    rescoring_features : dict
        A dictionary where keys are sources and values are lists of features.
    save_path : str, optional
        If provided, saves the plot to the specified path.
    sort : bool, optional
        If True, sorts the features by their importance in descending order.
        Default is False.
    error : bool, optional
        If True, adds error bars to the plot. Default is False.
    **kwargs : dict
        Additional plotting parameters:
        - 'figsize' : tuple, default (15, 10)
            Figure size in inches (width, height).
        - 'dpi' : int, default 300
            Resolution in dots per inch.
        - 'palette' : str, default 'crest'
            Seaborn color palette name. Options include 'crest', 'flare', 'mako',
            'rocket', 'tab10', 'husl', 'Set2', etc.

    Notes
    -----
    The function automatically detects the model type based on the presence of the corresponding attribute.
    For linear models, it uses hatch patterns to differentiate between positive and negative coefficients.
    For XGBoost models, it uses solid bars since the importances are always positive.

    The color palette is automatically scaled to match the number of feature sources, ensuring
    consistent colors between the bars and legend.

    Examples
    --------
    >>> # Use default crest palette
    >>> plot_feature_importance(models, rescoring_features, save_path='importance.png')

    >>> # Use a different palette
    >>> plot_feature_importance(models, rescoring_features, palette='flare', sort=True, error=True)
    """
    # Determine the model type based on the first model in the list.
    if hasattr(models[0].estimator, "coef_"):
        model_type = "linear"
    elif hasattr(models[0].estimator, "feature_importances_"):
        model_type = "xgb"
    else:
        raise ValueError(
            "Model type not recognized. Model must have 'estimator.coef_' for linear models or "
            "'estimator.feature_importances_' for XGBoost models."
        )

    if model_type == "linear":
        feature_importances = []
        for model in models:
            coefficients = model.estimator.coef_
            feature_importances.append(np.abs(coefficients).mean(axis=0))
            logger.debug(f"Model coefficients shape: {coefficients.shape}")

        average_feature_importance = np.mean(feature_importances, axis=0)
        std_feature_importance = np.std(feature_importances, axis=0)
        feature_signs = np.mean([model.estimator.coef_.mean(axis=0) for model in models], axis=0)

    elif model_type == "xgb":
        feature_importances = []
        for model in models:
            # Use the XGBoost feature importances directly as they are always positive
            imp = model.estimator.feature_importances_
            feature_importances.append(imp)
            logger.debug(f"Model feature importances shape: {imp.shape}")

        average_feature_importance = np.mean(feature_importances, axis=0)
        std_feature_importance = np.std(feature_importances, axis=0)
        feature_signs = np.ones_like(average_feature_importance)

    logger.debug(f"Total rescoring features: {len(sum(rescoring_features.values(), []))}")
    logger.debug(f"Average feature importance length: {len(average_feature_importance)}")
    logger.debug(f"Features: {sum(rescoring_features.values(), [])}")

    # Extract plotting parameters
    figsize = kwargs.get("figsize", (15, 10))
    dpi = kwargs.get("dpi", 300)
    palette_name = kwargs.get("palette", "Set2")

    all_features = []
    all_importances = []
    all_errors = []
    all_colors = []
    all_hatches = []  # Hatch patterns will be applied only for linear models.

    n_sources = len(rescoring_features)
    colors = sns.color_palette(palette_name, n_colors=n_sources)
    source_colors = dict(zip(rescoring_features.keys(), colors))

    for source, features in rescoring_features.items():
        color = source_colors[source]
        indices = [
            i for i, name in enumerate(sum(rescoring_features.values(), [])) if name in features
        ]
        source_importances = average_feature_importance[indices]
        source_std = std_feature_importance[indices]

        if model_type == "linear":
            source_signs = feature_signs[indices]

        if sort:
            sorted_indices = np.argsort(-source_importances)
        else:
            sorted_indices = np.arange(len(features))

        sorted_features = [features[i] for i in sorted_indices]
        sorted_importances = source_importances[sorted_indices]
        sorted_std = source_std[sorted_indices]

        all_features.extend(sorted_features)
        all_importances.extend(sorted_importances)
        all_errors.extend(sorted_std)
        all_colors.extend([color] * len(sorted_features))

        if model_type == "linear":
            # For linear models, use hatch patterns to differentiate positive and negative coefficients.
            # An empty hatch ('') for positive and '\\' for negative coefficients.
            sorted_signs = source_signs[sorted_indices]
            all_hatches.extend(["" if sign >= 0 else "\\\\" for sign in sorted_signs])
        else:
            all_hatches.extend([""] * len(sorted_features))

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    if error:
        bars = ax.barh(all_features, all_importances, xerr=all_errors, color=all_colors, capsize=5)
    else:
        bars = ax.barh(all_features, all_importances, color=all_colors)

    if model_type == "linear":
        for bar, hatch in zip(bars, all_hatches):
            bar.set_hatch(hatch)
        legend_hatches = [
            Patch(facecolor="white", edgecolor="black", hatch="", label="Positive"),
            Patch(facecolor="white", edgecolor="black", hatch="\\\\", label="Negative"),
        ]
        legend_colors = [
            Patch(facecolor=source_colors[source], edgecolor="black", label=source)
            for source in rescoring_features.keys()
        ]
        ax.legend(handles=legend_hatches + legend_colors, loc="best")
    else:
        legend_colors = [
            Patch(facecolor=source_colors[source], edgecolor="black", label=source)
            for source in rescoring_features.keys()
        ]
        ax.legend(handles=legend_colors, loc="best")

    ax.set_xlabel("Average Feature Importance")
    ax.set_ylabel("Feature")

    save_or_show_plot(save_path, logger)


def visualize_feature_correlation(psms: PsmContainer, save_path=None, **kwargs):
    """
    Visualize the correlation between features in a DataFrame using a scatter plot heatmap.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing the features to visualize.
    save_path : str, optional
        The file path to save the plot. If not provided, the plot is displayed.
    **kwargs : dict
        Additional plotting parameters such as `figsize`, `dpi`, `height`, etc.
    """
    rescoring_features = [item for sublist in psms.rescoring_features.values() for item in sublist]
    n_features = len(rescoring_features)

    default_height = max(8, min(20, 8 + n_features * 0.15))
    height = kwargs.get("height", default_height)
    dpi = kwargs.get("dpi", 300)

    corr = psms.psms[rescoring_features].corr()
    corr_mat = corr.stack().reset_index(name="correlation")

    g = sns.relplot(
        data=corr_mat,
        x="level_0",
        y="level_1",
        hue="correlation",
        size="correlation",
        palette="vlag",
        hue_norm=(-1, 1),
        edgecolor=".7",
        height=height,
        sizes=(50, 250),
        size_norm=(-0.2, 0.8),
        **{k: v for k, v in kwargs.items() if k not in ["height", "dpi", "figsize"]},
    )

    g.set(xlabel="", ylabel="", aspect="equal")
    g.despine(left=True, bottom=True)
    g.ax.margins(0.02)

    for label in g.ax.get_xticklabels():
        label.set_rotation(90)

    if n_features > 30:
        fontsize = max(6, 10 - n_features * 0.05)
        g.ax.tick_params(labelsize=fontsize)

    g.figure.suptitle("Feature Correlation Matrix", y=1.01, fontsize=14)
    plt.tight_layout()
    g.figure.dpi = dpi

    save_or_show_plot(save_path, logger)
