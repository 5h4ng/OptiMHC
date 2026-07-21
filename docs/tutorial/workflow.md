# Pipeline Workflow

This page explains the complete OptiMHC pipeline from input to output. Understanding the workflow helps you choose the right configuration for your data and troubleshoot issues.

## Overview

The pipeline executes in five sequential stages:

```
1. Configuration loading
2. Input parsing
3. Feature generation
4. Rescoring
5. Output & visualization
```

## Stage 1: Configuration

OptiMHC uses a layered configuration system. Settings are resolved in the following precedence order (highest to lowest):

1. **CLI flags** — command-line arguments override everything.
2. **YAML file** — values from the `--config` file.
3. **Default config** — built-in defaults for all settings.

The default configuration provides sensible starting values:

```yaml
outputDir: ./results
inputType: pepxml
decoyPrefix: DECOY_
visualization: true
saveModels: true
keepIntermediate: true
toFlashLFQ: true
numProcesses: 4
logLevel: INFO
rescore:
  testFDR: 0.01
  trainFDR: 0.01
  model: Percolator
  numJobs: 1
  seed: 1
```

The `Config` class deep-merges your YAML file with these defaults, so you only need to specify what differs.

## Stage 2: Input Parsing

OptiMHC accepts two input formats:

### PepXML

The PepXML parser extracts PSMs from the XML structure produced by search engines (e.g., Comet, X!Tandem). For each PSM it extracts:

- Spectrum metadata (scan number, spectrum ID, charge, retention time)
- Mass values (experimental and calculated neutral mass)
- Peptide sequence with modifications
- Protein accessions
- All search engine scores

The parser preserves every `search_hit`, including its candidate rank, then computes the existing derived search features. See [Original Features](features/original.md) for details.

### PIN (Percolator Input)

The PIN parser reads one tab-separated Percolator input file per call. Columns that are not recognized metadata are declared directly as rescoring features.

### PsmContainer

Regardless of input format, parsing produces a `PsmContainer`: a small wrapper around one canonical pandas DataFrame plus an immutable tuple of feature column names.

```python
psms.df[["run", "scan", "rank", "sequence", "charge"]]
psms.feature_columns
```

The canonical identity of a physical spectrum is `run` plus `scan`; `rank` identifies competing candidates and is always a declared model feature. Multi-file combination is handled by the pipeline in the configured file order.

## Stage 3: Feature Generation

The pipeline iterates over the `featureGenerator` list from the configuration. For each entry, it:

1. Instantiates the feature class by name.
2. Calls `generate_features()`, which returns a DataFrame.
3. Attaches explicitly named numeric columns with `add_features(..., on=..., columns=...)` using a canonical key such as `psm_id`, `run` plus `scan`, or `sequence`.

After all generators have run, the PsmContainer holds the complete feature matrix. The
orchestration layer also records a run-local manifest mapping `Original` and each generator
name to its declared output columns. This grouping is not stored as mutable state in the
PsmContainer.

Available generators are documented in detail in the [Features](features/index.md) section:

| Feature            | Typical canonical join key |
| ------------------ | -------------------------- |
| Basic              | `psm_id`                   |
| SpectralSimilarity | `psm_id`                   |
| DeepLC             | `psm_id`                   |
| OverlappingPeptide | `sequence`                 |
| PWM                | `sequence`                 |
| MHCflurry          | `sequence`                 |
| NetMHCpan          | `sequence`                 |
| NetMHCIIpan        | `sequence`                 |

## Stage 4: Rescoring

Rescoring uses the [mokapot](https://mokapot.readthedocs.io/) framework. The pipeline:

1. **Builds a mokapot dataset** — the same PIN projection is used for in-process rescoring and PIN serialization. Mokapot groups candidates by `filename` (`run`) plus `ScanNr` (`scan`), while rank remains a feature.
2. **Trains a model** — `mokapot.brew()` performs semi-supervised learning with 3-fold cross-validation:
   - Trains the model on the training fold.
   - Scores PSMs in the test fold.
   - Repeats for all folds.
3. **Assigns q-values** using target-decoy competition at the specified `testFDR`. The configured integer seed (default `1`) is passed to the model and `mokapot.brew`.

### Available Models

| Model            | Description                                                                                                                                                                       |
| ---------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Percolator**   | Linear SVM (default). Fast and robust. Uses `mokapot.PercolatorModel`.                                                                                                            |
| **XGBoost**      | Gradient-boosted trees. Hyperparameters tuned via `GridSearchCV` (3-fold CV, ROC-AUC). Searches over `scale_pos_weight`, `max_depth`, `min_child_weight`, `gamma`.                |
| **RandomForest** | Random forest classifier. Hyperparameters tuned via `GridSearchCV` (3-fold CV, ROC-AUC). Searches over `class_weight`, `max_depth`, `min_samples_split`, `min_impurity_decrease`. |

## Stage 5: Output & Visualization

### Output Files

- **mokapot result files** — PSM-level and peptide-level results with scores and q-values.
- **PIN file** — the complete feature matrix in Percolator input format (useful for downstream tools).
- **Models** — serialized rescoring models (when `saveModels: true`).
- **FlashLFQ file** — accepted peptides in quantification-ready format (when
  `toFlashLFQ: true`). Retention time is written in minutes; missing calculated masses are
  derived from the canonical peptidoform.
- **BA Parquet intermediate results** — best binding prediction per peptide and predictor, documented under [Binding-Affinity Intermediate Results](features/binding-affinity.md#binding-affinity-intermediate-results).

### Visualizations

When `visualization: true`, the pipeline produces:

| Plot                           | Description                                                                                              |
| ------------------------------ | -------------------------------------------------------------------------------------------------------- |
| **qvalues.png**                | Number of PSMs and peptides accepted at each q-value threshold.                                          |
| **feature_importance.png**     | Bar chart showing the weight or importance of each feature in the trained model.                         |
| **feature_correlation.png**    | Heatmap of pairwise Pearson correlations among all rescoring features.                                   |
| **target_decoy_histogram.png** | KDE histograms comparing the distribution of each feature for targets vs. decoys (top-ranked hits only). |

## Experiment Mode

Experiment mode (`optimhc experiment --config ...`) shares stages 1–3 with the standard pipeline but then runs multiple rescoring experiments in parallel. Each experiment normally selects complete generator categories through `source`; advanced experiments may instead select exact DataFrame columns through `features`. An experiment cannot use both.

This is useful for:

- **Ablation studies** — measuring the contribution of each feature group.
- **Model comparison** — comparing Percolator vs. XGBoost vs. RandomForest on the same data.

Each experiment runs in its own process and writes results to a separate subdirectory. A shared PIN file and feature correlation plot are generated once for the full feature set.
