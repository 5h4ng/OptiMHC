# Features Overview

OptiMHC uses a modular feature generation system. Each feature computes a set of scoring columns that are added to the PSM data before rescoring. By combining features from multiple orthogonal sources, the rescoring model can better distinguish true peptide identifications from false ones.

## Feature Architecture

All features inherit from `BaseFeatureGenerator` and implement a common interface:

- **`feature_columns`** — the list of output column names.
- **`feature_groups(name)`** — the experiment category and columns exposed by the generator.
- **`id_column`** — the key column(s) used to merge features back into the PsmContainer.
- **`generate_features()`** — computes and returns a DataFrame of features.

Features are configured in the YAML file as a list of `{name, params}` entries. The pipeline instantiates each feature, calls `generate_features()`, and merges the result into the PsmContainer.
The same declaration builds a run-local feature manifest, allowing experiment YAML to select a
generator name without listing all of its output columns.

## Feature Categories

| Category | Name | Columns | Description |
|---|---|---|---|
| [Original Features](original.md) | PepXML / PIN parser | Variable | Search engine scores and derived mass/ion features |
| [Basic Features](basic.md) | `Basic` | 5 | Peptide sequence properties (length, entropy, AA composition) |
| [Spectral Similarity](spectral-similarity.md) | `SpectralSimilarity` | 8 | Similarity between experimental and predicted spectra |
| [Retention Time Deviation](rt-deviation.md) | `DeepLC` | 5 | Deviation between observed and predicted retention time |
| [Binding Affinity](binding-affinity.md) | `NetMHCpan`, `NetMHCIIpan`, `MHCflurry` | 3–4 per tool | MHC binding affinity predictions |
| [PWM Score](pwm.md) | `PWM` | 1–3 per allele | Position weight matrix binding score |
| [Overlapping Peptide Score](overlapping.md) | `OverlappingPeptide` | 4 | Graph-based peptide overlap and contig assembly features |

## Peptide representation

PIN and pepXML readers provide the unmodified peptide sequence to feature
generators. Modifications and their positions are stored separately for features
that need them.

??? note "Selenocysteine (U) handling"
    Selenocysteine (`U`) is an extremely rare amino acid that most prediction tools do not support. During preprocessing, `U` is automatically replaced with cysteine (`C`) to ensure compatibility with external tools such as Koina, DeepLC, MHCflurry, NetMHCpan, and the PWM scoring matrices.

## Missing Value Handling

When a feature cannot be computed for a PSM (e.g., peptide length outside the supported range for a binding predictor), the missing values are filled with the **median** of the successfully computed values for that column. This ensures no NaN values propagate to the rescoring model.
