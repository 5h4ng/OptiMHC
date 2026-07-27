import logging

import matplotlib.pyplot as plt
import seaborn as sns

from optimhc.psm_container import PsmContainer
from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)


def visualize_target_decoy_features(psms: PsmContainer, num_cols=5, save_path=None, **kwargs):
    """
    Visualize the distribution of features in a DataFrame using kernel density estimation plots.

    Parameters
    ----------
    psms : PsmContainer
        A PsmContainer object containing the features to visualize.
    num_cols : int, optional
        The number of columns in the plot grid. Default is 5.
    save_path : str, optional
        The file path to save the plot. If not provided, the plot is displayed.
    **kwargs : dict
        Additional plotting parameters such as `figsize` and `dpi`, etc.

    Notes
    -----
    This function:
    1. Extracts rescoring features from the PsmContainer
    2. Filters out features with only one unique value
    3. Creates a grid of plots showing the distribution of each feature
    4. Separates target and decoy PSMs in each plot
    5. Uses kernel density estimation to show the distribution shape
    """
    rescoring_features = [feature for feature in psms.feature_columns if feature != "rank"]

    # drop features that only have one value
    rescoring_features = [
        feature for feature in rescoring_features if psms.df[feature].nunique() > 1
    ]

    num_features = len(rescoring_features)
    num_rows = (num_features + num_cols - 1) // num_cols

    figsize = kwargs.get("figsize", (15, num_rows * 15 / num_cols))
    dpi = kwargs.get("dpi", 300)

    fig, axes = plt.subplots(num_rows, num_cols, figsize=figsize, dpi=dpi)
    axes = axes.flatten()

    psms_top_hits = psms.df[psms.df["rank"] == 1].copy()
    num_true_hits = len(psms_top_hits[~psms_top_hits["is_decoy"]])
    num_decoys = len(psms_top_hits[psms_top_hits["is_decoy"]])
    logger.debug(f"Number of true hits: {num_true_hits}")
    logger.debug(f"Number of decoys: {num_decoys}")
    psms_top_hits["target_decoy"] = psms_top_hits["is_decoy"].map({False: "Target", True: "Decoy"})

    for i, feature in enumerate(rescoring_features):
        try:
            ax = axes[i]
            sns.histplot(
                data=psms_top_hits,
                x=feature,
                hue="target_decoy",
                ax=ax,
                bins="auto",
                common_bins=True,
                multiple="dodge",
                fill=True,
                alpha=0.3,
                stat="frequency",
                kde=True,
                linewidth=0,
            )
            ax.set_title(feature)
            ax.set_xlabel("")
            ax.set_ylabel("")

        except Exception as e:
            logger.error(f"Error plotting feature {feature}: {e}")
            ax.set_visible(False)

    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    save_or_show_plot(save_path, logger)
