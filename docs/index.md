# OptiMHC

**An optimized rescoring pipeline for immunopeptidomics data that significantly enhances peptide identification performance.**

OptiMHC integrates multiple rescoring features with machine learning-based rescoring to maximize the number of confidently identified peptides from mass spectrometry experiments.

## How It Works

```
Input files (PepXML / PIN)
  → Parsing & feature extraction
    → PsmContainer (central data structure)
      → Feature generation (Basic, Spectral, RT, MHC binding, PWM, Overlap, …)
        → Machine learning rescoring (Percolator / XGBoost / RandomForest)
          → Visualization & output
```

1. **Parse** search engine results from PepXML or PIN format into a unified `PsmContainer`.
2. **Generate features** using a configurable set of features — each adds new scoring dimensions to the PSM data.
3. **Rescore** PSMs with machine learning models (Percolator SVM, XGBoost, or RandomForest) trained via mokapot to separate targets from decoys at a controlled FDR.
4. **Visualize** results with q-value curves, feature importance plots, and target/decoy distributions.

## Getting Started

- [Installation](getting-started/installation.md) — set up OptiMHC on your system
- [Quick Start](getting-started/quickstart.md) — run your first rescoring pipeline in minutes

## Learn More

- [Tutorial](tutorial/index.md) — examples, pipeline walkthrough, and feature explanations
- [API Reference](api/index.md) — detailed module and class documentation
- [Development](development/index.md) — set up a development environment
