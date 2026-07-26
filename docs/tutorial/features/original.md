# Original Features

Original features are the initial scores and derived quantities declared by an input reader.
They are ordinary PSM DataFrame columns listed in `psms.feature_columns`. Experiment mode
exposes them as the `Original` entry in its run-local feature manifest; the PsmContainer itself
does not maintain a source registry. The exact columns depend on the input format.

## PepXML Features

When parsing PepXML files, OptiMHC extracts raw search engine scores and computes several derived features.

### Search Engine Score

All `<search_score>` elements from the PepXML file are extracted as features. The exact columns depend on the search engine that produced the file. For example, Comet produces:

- `xcorr` — cross-correlation score
- `deltacn` — delta correlation (difference to next-best match)
- `deltacnstar` — delta CN star
- `spscore` — preliminary score
- `sprank` — SP rank
- `expect` — expectation value (E-value)

Other search engines (X!Tandem, MSFragger, etc.) produce their own score columns, all of which are automatically included.

### Mass Difference Features

Let \( m_\mathrm{exp} \) and \( m_\mathrm{calc} \) denote the experimental and calculated precursor neutral masses, respectively. The mass difference \( \Delta m \) is:

\[
\Delta m = m_\mathrm{exp} - m_\mathrm{calc}
\]

The experimental and calculated m/z values are derived from the charge state \( z \) using the proton mass \( m_\mathrm{H} = 1.00727646677 \) Da:

\[
(m/z)_\mathrm{exp} = \frac{m_\mathrm{exp}}{z} + m_\mathrm{H}, \quad
(m/z)_\mathrm{calc} = \frac{m_\mathrm{calc}}{z} + m_\mathrm{H}
\]

The absolute m/z difference \( \Delta_{mz} \) is:

\[
\Delta_{mz} = \left|(m/z)_\mathrm{exp} - (m/z)_\mathrm{calc}\right|
\]

### Matched Ion Ratio

When the search engine reports matched and total ions, the matched ion ratio \( R_\mathrm{ion} \) is computed:

\[
R_\mathrm{ion} = \frac{n_\mathrm{matched}}{n_\mathrm{total}}
\]

This feature is only included if both values are present and \( n_\mathrm{total} > 0 \) for all PSMs.

### Charge One-Hot Encoding

The precursor charge state is one-hot encoded into binary columns:

- `charge_1`, `charge_2`, `charge_3`, ...

For a PSM with charge 2, `charge_2 = 1` and all other charge columns are 0. This allows the rescoring model to learn charge-specific effects without treating charge as a continuous variable.

### Log Transformation of P-values and E-values

Many search engine scores span several orders of magnitude (e.g., E-values from \( 10^{-20} \) to \( 10^{2} \)). Following the approach used in mokapot, OptiMHC applies automatic log-transformation to compress these ranges into a scale more suitable for machine learning models.

**Scientific notation values:** If values contain scientific notation and span 4 or more orders of magnitude, the log transform for a value \( x = a \times 10^{b} \) is:

\[
\log_{10}(x) = \log_{10}(a) + b
\]

For zero values, the log is set to one less than the minimum observed log value.

**Numeric columns:** For non-negative, non-binary numeric columns where \( x_\mathrm{max} / x_\mathrm{min}^{+} \geq 10{,}000 \) (where \( x_\mathrm{min}^{+} \) is the smallest nonzero value):

\[
x' = \begin{cases}
\log_{10}(x) & \text{if } x > 0 \\
\min\!\big(\log_{10}(x_\mathrm{min}^{+})\big) - 1 & \text{if } x = 0
\end{cases}
\]

The `num_matched_peptides` column is always log-transformed as \( \log_{10}(x) \).

### Decoy Detection

A PSM is labeled as a **decoy** if the first protein accession starts with the `decoyPrefix` (default: `DECOY_`). If any alternative protein accession does *not* start with the prefix, the PSM is re-labeled as a **target**.

### Additional Metadata

The reader also stores these PSM fields. Only `rank` is used as a rescoring feature:

- `run` — acquisition run / raw-file identifier
- `scan` — scan number
- `rank` — candidate rank and rescoring feature
- `is_decoy` — Boolean target/decoy state
- `calc_mass` — calculated neutral mass
- `sequence`, `mods`, `mod_sites` — normalized peptidoform fields
- `proteins` — protein accessions
- `charge` — internal precursor charge metadata; input charge columns remain features
- `retention_time` — seconds for pepXML; PIN values are multiplied by 60 when the column maximum is below 500

---

## PIN Features

When parsing PIN files, the feature set is determined by the file itself. Columns not recognized
as metadata are retained as rescoring features when all their values are numeric. Non-numeric
extra columns are skipped with a warning. Rank is preserved from a rank column when present,
otherwise parsed from a trailing `_<rank>` in `SpecId`, and synthesized as `1` only when neither
representation is available. See [PIN Input](../pin-format.md) for the supported metadata,
peptide, and charge formats.

Column name matching is case-insensitive. Charge columns matching the pattern `charge[_]?\d+` (e.g., `charge1`, `charge_2`) are detected automatically. For each PSM, the charge value is determined by which charge column contains a value of 1.

The Label column uses the convention `1` for target and `-1` for decoy.
