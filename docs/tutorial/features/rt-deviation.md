# Retention Time Deviation

The **DeepLC** feature (`name: DeepLC`) predicts peptide retention times using the [DeepLC](https://github.com/compomics/DeepLC) deep learning model and computes the deviation between predicted and observed retention times. Correct peptide identifications should show smaller deviations, while incorrect matches tend to have larger discrepancies.

**Source name:** `DeepLC`

## Output Columns

| Column | Description |
|---|---|
| `observed_retention_time` | Retention time recorded in the input data (seconds) |
| `predicted_retention_time` | Retention time predicted by DeepLC |
| `retention_time_diff` | Signed difference: predicted minus observed |
| `abs_retention_time_diff` | Absolute value of the retention time difference |
| `retention_time_ratio` | Ratio of the smaller to the larger of predicted and observed |

## Computation

### Step 1: Preprocessing

Peptide sequences are preprocessed for DeepLC input:

1. **Flanking amino acids** are stripped (e.g., `K.PEPTIDE.R` → `PEPTIDE`).
2. **Modifications** are converted from mass-annotated format to UNIMOD notation using the `modificationMap`. For example, `M[147.035]` is converted to a UNIMOD-indexed modification string that DeepLC understands.

### Step 2: Calibration

DeepLC benefits from calibration on high-confidence PSMs from the same LC-MS run. The calibration procedure:

1. **Sort** all PSMs by the `calibrationCriteria` column (a search engine score). If `lowerIsBetter` is true, sort ascending; otherwise, sort descending.
2. **Select** the top PSMs as the calibration set:
    - If `calibrationSize` is a float (e.g., 0.1), take the top 10% of PSMs.
    - If `calibrationSize` is an integer (e.g., 500), take the top 500 PSMs.
3. **Filter** to target PSMs only (decoys are excluded from calibration).
4. **Calibrate** the DeepLC predictor using the observed retention times of the calibration set.

### Step 3: Prediction

The calibrated DeepLC model predicts retention times for all PSMs in the dataset.

### Step 4: Feature Computation

Let \( t_\mathrm{obs} \) and \( t_\mathrm{pred} \) denote the observed and predicted retention times, respectively.

The signed difference \( \Delta t \), its absolute value \( |\Delta t| \), and the retention time ratio \( R_t \) are:

\[
\Delta t = t_\mathrm{pred} - t_\mathrm{obs}
\]

\[
|\Delta t| = |t_\mathrm{pred} - t_\mathrm{obs}|
\]

\[
R_t = \frac{\min(t_\mathrm{pred},\; t_\mathrm{obs})}{\max(t_\mathrm{pred},\; t_\mathrm{obs})}
\]

\( R_t \) is bounded in \((0, 1]\) and equals 1 when the predicted and observed times are identical.

### Missing Values

Any PSMs for which DeepLC cannot produce a prediction are filled with the **median** of each feature column.

## Configuration

```yaml
featureGenerator:
  - name: DeepLC
    params:
      calibrationCriteria: expect    # Column name used for selecting calibration PSMs
      lowerIsBetter: true            # Whether lower values of the criteria are better
      calibrationSize: 0.1           # Fraction (float) or count (int) of top PSMs for calibration
```

| Parameter | Default | Description |
|---|---|---|
| `calibrationCriteria` | *(required)* | Name of a search engine score column to rank PSMs for calibration |
| `lowerIsBetter` | `false` | Set to `true` if lower values of the calibration criteria indicate better PSMs (e.g., E-values) |
| `calibrationSize` | `0.15` | Fraction of PSMs (float) or absolute count (int) for the calibration set |

