import logging

import matplotlib.pyplot as plt

from optimhc.visualization.save_or_show_plot import save_or_show_plot

logger = logging.getLogger(__name__)


def plot_qvalues(
    results,
    save_path=None,
    dpi=300,
    figsize=(15, 10),
    threshold=0.05,
    colors=None,
    **kwargs,
):
    """
    Plot q-values for the given results.

    Parameters
    ----------
    results : object or list
        A list of results objects or a single result object.
        Each result object should have a method `plot_qvalues`.
    save_path : str, optional
        If provided, saves the plot to the specified path.
    dpi : int, optional
        The resolution of the plot. Default is 300.
    figsize : tuple, optional
        The size of the figure. Default is (15, 10).
    threshold : float, optional
        The q-value threshold for plotting. Default is 0.05.
    colors : list, optional
        A list of colors for the plots. If not provided, uses default colors.
    **kwargs : dict
        Additional plotting parameters.

    Returns
    -------
    None
        The function displays or saves the plot.

    Notes
    -----
    This function:
    1. Creates a figure with two subplots for PSMs and peptides
    2. Plots q-values for each result with different colors
    3. Adds legends and titles to each subplot
    4. Saves or displays the plot based on save_path
    """
    if not isinstance(results, list):
        results = [results]

    if colors is None:
        colors = [
            "#1f77b4",
            "#ff7f0e",
            "#2ca02c",
            "#d62728",
            "#9467bd",
            "#8c564b",
            "#e377c2",
            "#7f7f7f",
            "#bcbd22",
            "#17becf",
        ]

    fig, axs = plt.subplots(1, 2, figsize=figsize, dpi=dpi)

    for i, result in enumerate(results):
        for ax, level in zip(axs, ["psms", "peptides"]):
            result.plot_qvalues(
                level=level,
                c=colors[i % len(colors)],
                ax=ax,
                threshold=threshold,
                label=f"Result {i + 1}" if len(results) > 1 else "Results",
                linewidth=1,
                **kwargs,
            )
            ax.legend(frameon=False)
            ax.set_title(f"{level}")

    plt.tight_layout()
    return save_or_show_plot(save_path, logger)
