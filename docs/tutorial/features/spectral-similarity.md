# Spectral Similarity

The **SpectralSimilarity** feature (`name: SpectralSimilarity`) compares experimentally observed MS2 spectra against theoretical spectra predicted by a deep learning model served via [Koina](https://koina.wilhelmlab.org/). The similarity between observed and predicted fragmentation patterns is a powerful indicator of correct peptide-spectrum matches.

**Source name:** `SpectralSimilarity`

## Output Columns

| Column | Description |
|---|---|
| `spectral_angle_similarity` | Normalized spectral angle between experimental and predicted spectra |
| `cosine_similarity` | Cosine similarity between intensity vectors |
| `pearson_correlation` | Pearson correlation coefficient |
| `spearman_correlation` | Spearman rank correlation coefficient |
| `mean_squared_error` | Mean squared error on L2-normalized intensity vectors |
| `unweighted_entropy_similarity` | Entropy-based similarity measure |
| `predicted_seen_nonzero` | Count of predicted peaks that were observed in the experimental spectrum |
| `predicted_not_seen` | Count of predicted peaks that were not observed |

## Computation

### Step 1: Spectrum Extraction

Experimental MS2 spectra are extracted from mzML files. Each spectrum is associated with a PSM through the spectrum ID pattern. The m/z and intensity arrays are sorted by m/z.

### Step 2: Theoretical Spectrum Prediction

Peptide sequences (with modifications mapped to UNIMOD notation via `modificationMap`) and charge states are sent to a Koina server, which returns predicted fragment ion m/z values and intensities using the specified deep learning model (e.g., `AlphaPeptDeep_ms2_generic`).

### Step 3: Peak Alignment

For each predicted peak at m/z value \( m \), a tolerance window is computed in parts-per-million (ppm):

\[
m_\mathrm{low} = m \cdot \left(1 - \frac{\tau}{10^6}\right), \quad
m_\mathrm{high} = m \cdot \left(1 + \frac{\tau}{10^6}\right)
\]

where \( \tau \) is the ppm tolerance (default 20). Within this window, the experimental peak with the **highest intensity** is selected as the match. If no experimental peak falls within the window, the aligned experimental intensity is set to 0.

This produces an aligned experimental intensity vector \( \mathbf{e} \) of the same length as the predicted intensity vector \( \mathbf{p} \).

### Step 4: Top-N Peak Selection

Only the top \( N \) peaks by predicted intensity are retained for similarity computation (default \( N = 36 \)). This focuses the comparison on the most informative fragment ions.

### Step 5: Similarity Metrics

All metrics are computed on the aligned, top-N intensity vectors \( \mathbf{e} \) (experimental) and \( \mathbf{p} \) (predicted).

#### Spectral Angle Similarity

The spectral angle similarity \( S_\mathrm{SA} \) is derived from the cosine of the angle \( \theta \) between the two vectors:

\[
\cos\theta = \frac{\mathbf{e} \cdot \mathbf{p}}{\|\mathbf{e}\| \, \|\mathbf{p}\|}
\]

\[
S_\mathrm{SA} = 1 - \frac{2 \arccos(\cos\theta)}{\pi}
\]

\( S_\mathrm{SA} \) ranges from 0 (orthogonal spectra) to 1 (identical spectra). The \( 2/\pi \) normalization maps the angle from \([0, \pi/2]\) to \([0, 1]\).

#### Cosine Similarity

The cosine similarity \( S_{\mathrm{cos}} \) is:

\[
S_{\mathrm{cos}} = \frac{\mathbf{e} \cdot \mathbf{p}}{\|\mathbf{e}\| \, \|\mathbf{p}\|}
\]

#### Pearson Correlation

The Pearson correlation coefficient \( r \) is:

\[
r = \frac{\sum_{i} (e_i - \bar{e})(p_i - \bar{p})}{\sqrt{\sum_{i} (e_i - \bar{e})^2} \sqrt{\sum_{i} (p_i - \bar{p})^2}}
\]

where \( \bar{e} \) and \( \bar{p} \) are the means of the experimental and predicted vectors.

#### Spearman Correlation

The Spearman correlation \( \rho \) is the Pearson correlation computed on the **ranks** of the values (using average ranking for ties) rather than the values themselves.

#### Mean Squared Error

The MSE is computed on L2-normalized vectors. Let \( \hat{e}_i = e_i / \|\mathbf{e}\| \) and \( \hat{p}_i = p_i / \|\mathbf{p}\| \). Then:

\[
\mathrm{MSE} = \frac{1}{n} \sum_{i=1}^{n} (\hat{e}_i - \hat{p}_i)^2
\]

#### Unweighted Entropy Similarity

First, normalize intensities to probability distributions \( \mathbf{q}^{(e)} \) and \( \mathbf{q}^{(p)} \):

\[
q_i^{(e)} = \frac{e_i}{\sum_j e_j}, \quad q_i^{(p)} = \frac{p_i}{\sum_j p_j}
\]

Compute the Shannon entropy of each distribution and their mixture \( \mathbf{m} \):

\[
S_e = -\sum_i q_i^{(e)} \ln q_i^{(e)}, \quad
S_p = -\sum_i q_i^{(p)} \ln q_i^{(p)}
\]

\[
m_i = \frac{1}{2}\left(q_i^{(e)} + q_i^{(p)}\right), \quad
S_m = -\sum_i m_i \ln m_i
\]

The entropy similarity \( S_\mathrm{ent} \) is:

\[
S_\mathrm{ent} = 1 - \frac{2 S_m - S_e - S_p}{\ln 4}
\]

The numerator is related to the Jensen-Shannon divergence. Dividing by \( \ln 4 \) normalizes the result to \([0, 1]\).

#### Peak Matching Counts

- **`predicted_seen_nonzero`** — the number of predicted peaks (with intensity > 0) that have a matched experimental peak with intensity > 0.
- **`predicted_not_seen`** — the number of predicted peaks (with intensity > 0) that have no matching experimental peak (aligned intensity = 0).

## Configuration

```yaml
featureGenerator:
  - name: SpectralSimilarity
    params:
      mzmlDir: ../data               # Directory containing mzML files
      spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+  # Regex to link PSMs to mzML files
      model: AlphaPeptDeep_ms2_generic          # Koina prediction model
      collisionEnergy: 28            # Collision energy for prediction
      instrument: LUMOS              # Instrument type for prediction
      tolerance: 20                  # Peak matching tolerance in ppm
      numTopPeaks: 36                # Number of top peaks to compare
      url: koina.wilhelmlab.org:443  # Koina server gRPC endpoint
      ssl: true                      # Use SSL for gRPC connection
```

!!! note
    For a local Koina/Triton server, set `url: 127.0.0.1:8500` and `ssl: false`. The default public endpoint is `koina.wilhelmlab.org:443` with SSL enabled.

## Requirements

This feature requires **mzML files** (raw spectra) and access to a **Koina server** — either the [public endpoint](https://koina.wilhelmlab.org/) or a self-hosted [Triton Inference Server](https://github.com/wilhelm-lab/koina).
