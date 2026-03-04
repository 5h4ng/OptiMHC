import numba as nb
import numpy as np


@nb.njit(cache=True)
def align_peaks(
    exp_mz_sorted: np.ndarray,
    exp_intensity_sorted: np.ndarray,
    pred_mz_sorted: np.ndarray,
    tolerance_ppm: float,
):
    """
    Align sorted experimental peaks to sorted predicted peaks using ppm tolerance.

    For each predicted peak, find the experimental peak within the tolerance window
    that has the highest intensity.

    Returns
    -------
    aligned_exp_intensity : float64 array of length n_pred
    matched_exp_indices   : int64 array of length n_pred (-1 = no match)
    """
    n_pred = len(pred_mz_sorted)
    n_exp = len(exp_mz_sorted)

    aligned_exp_intensity = np.zeros(n_pred, dtype=np.float64)
    matched_exp_indices = np.full(n_pred, -1, dtype=np.int64)

    start_pos = 0

    for i in range(n_pred):
        pred_peak_mz = pred_mz_sorted[i]
        fragment_min = pred_peak_mz * (1.0 - tolerance_ppm / 1e6)
        fragment_max = pred_peak_mz * (1.0 + tolerance_ppm / 1e6)

        matched_int = 0.0
        past_start = 0

        while start_pos + past_start < n_exp:
            exp_peak_mz = exp_mz_sorted[start_pos + past_start]
            if exp_peak_mz < fragment_min:
                start_pos += 1
            elif exp_peak_mz <= fragment_max:
                exp_peak_int = exp_intensity_sorted[start_pos + past_start]
                if exp_peak_int > matched_int:
                    matched_int = exp_peak_int
                    matched_exp_indices[i] = start_pos + past_start
                past_start += 1
            else:
                break

        aligned_exp_intensity[i] = matched_int

    return aligned_exp_intensity, matched_exp_indices
