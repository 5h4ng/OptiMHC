import os

import matplotlib.pyplot as plt


def save_or_show_plot(save_path, logger, tight_layout=True):
    if tight_layout:
        plt.tight_layout()
    if save_path:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, bbox_inches="tight")
        logger.info(f"Plot saved to {save_path}")
    else:
        plt.show()
