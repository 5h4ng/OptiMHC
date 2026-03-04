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


@nb.njit(cache=True)
def _rankdata_average(arr):
    n = len(arr)
    sorter = np.argsort(arr)
    ranks = np.empty(n, dtype=np.float64)
    sorted_arr = arr[sorter]
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_arr[j] == sorted_arr[i]:
            j += 1
        avg = 0.5 * (i + 1 + j)
        for k in range(i, j):
            ranks[sorter[k]] = avg
        i = j
    return ranks


@nb.njit(cache=True)
def compute_similarity_features(exp_vector, pred_vector):
    """
    Compute all similarity metrics between two aligned intensity vectors.

    Returns
    -------
    tuple of 8 float64 values:
        (spectral_angle_similarity, cosine_similarity,
         pearson_correlation, spearman_correlation,
         mean_squared_error, unweighted_entropy_similarity,
         predicted_seen_nonzero, predicted_not_seen)
    """
    n = len(exp_vector)

    exp_sq = 0.0
    pred_sq = 0.0
    dot = 0.0
    exp_sum = 0.0
    pred_sum = 0.0
    for i in range(n):
        e = exp_vector[i]
        p = pred_vector[i]
        exp_sq += e * e
        pred_sq += p * p
        dot += e * p
        exp_sum += e
        pred_sum += p

    exp_l2 = np.sqrt(exp_sq)
    pred_l2 = np.sqrt(pred_sq)

    # spectral angle similarity
    if exp_l2 > 0.0 and pred_l2 > 0.0:
        cos_val = dot / (exp_l2 * pred_l2)
        cos_val = max(-1.0, min(1.0, cos_val))
        spectral_angle = 1.0 - 2.0 * np.arccos(cos_val) / np.pi
    else:
        spectral_angle = 0.0

    # cosine similarity
    if exp_sum == 0.0 or pred_sum == 0.0:
        cosine_sim = 0.0
    elif exp_l2 > 0.0 and pred_l2 > 0.0:
        cosine_sim = dot / (exp_l2 * pred_l2)
    else:
        cosine_sim = 0.0

    # MSE on L2-normalized vectors
    mse = 0.0
    for i in range(n):
        e = exp_vector[i] / exp_l2 if exp_l2 > 0.0 else 0.0
        p = pred_vector[i] / pred_l2 if pred_l2 > 0.0 else 0.0
        diff = e - p
        mse += diff * diff
    mse /= n if n > 0 else 1.0

    # pearson correlation
    exp_mean = exp_sum / n if n > 0 else 0.0
    pred_mean = pred_sum / n if n > 0 else 0.0
    num_p = 0.0
    exp_var = 0.0
    pred_var = 0.0
    for i in range(n):
        de = exp_vector[i] - exp_mean
        dp = pred_vector[i] - pred_mean
        num_p += de * dp
        exp_var += de * de
        pred_var += dp * dp
    if exp_var > 0.0 and pred_var > 0.0:
        pearson = num_p / np.sqrt(exp_var * pred_var)
    else:
        pearson = 0.0

    # spearman correlation
    if exp_var > 0.0 and pred_var > 0.0:
        exp_ranks = _rankdata_average(exp_vector)
        pred_ranks = _rankdata_average(pred_vector)
        rmean_e = 0.0
        rmean_p = 0.0
        for i in range(n):
            rmean_e += exp_ranks[i]
            rmean_p += pred_ranks[i]
        rmean_e /= n
        rmean_p /= n
        snum = 0.0
        svar_e = 0.0
        svar_p = 0.0
        for i in range(n):
            de = exp_ranks[i] - rmean_e
            dp = pred_ranks[i] - rmean_p
            snum += de * dp
            svar_e += de * de
            svar_p += dp * dp
        if svar_e > 0.0 and svar_p > 0.0:
            spearman = snum / np.sqrt(svar_e * svar_p)
        else:
            spearman = 0.0
    else:
        spearman = 0.0

    # unweighted entropy similarity
    s_exp = 0.0
    s_pred = 0.0
    s_mixed = 0.0
    for i in range(n):
        ep = exp_vector[i] / exp_sum if exp_sum > 0.0 else 0.0
        pp = pred_vector[i] / pred_sum if pred_sum > 0.0 else 0.0
        if ep > 0.0:
            s_exp -= ep * np.log(ep)
        if pp > 0.0:
            s_pred -= pp * np.log(pp)
        mp = 0.5 * (ep + pp)
        if mp > 0.0:
            s_mixed -= mp * np.log(mp)
    entropy_sim = 1.0 - (2.0 * s_mixed - s_exp - s_pred) / np.log(4.0)

    # predicted seen/not seen
    seen = 0.0
    not_seen = 0.0
    for i in range(n):
        if pred_vector[i] > 0.0:
            if exp_vector[i] > 0.0:
                seen += 1.0
            else:
                not_seen += 1.0

    return (
        spectral_angle,
        cosine_sim,
        pearson,
        spearman,
        mse,
        entropy_sim,
        seen,
        not_seen,
    )
