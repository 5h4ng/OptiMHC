# Quick Start

This guide walks you through running your first OptiMHC rescoring pipeline.

## 1. Prepare a Configuration File

Create a YAML file (e.g., `my_config.yaml`) that describes your input data and desired features:

```yaml
experimentName: my_first_run
inputType: pepxml
inputFile:
  - path/to/your/search_results.pep.xml
decoyPrefix: DECOY_
outputDir: ./results

allele:
  - HLA-A*02:01

featureGenerator:
  - name: Basic
  - name: OverlappingPeptide
    params:
      minOverlapLength: 7
      minLength: 7
      maxLength: 20
  - name: PWM
    params:
      class: I

rescore:
  testFDR: 0.01
  model: Percolator
```

!!! tip
    Start with lightweight features like **Basic**, **OverlappingPeptide**, and **PWM** — they require no external setup. You can add more features (SpectralSimilarity, DeepLC, MHCflurry, NetMHCpan) later.

## 2. Run the Pipeline

```bash
optimhc pipeline --config my_config.yaml
```

Or pass options directly on the command line:

```bash
optimhc pipeline \
  --inputType pepxml \
  --inputFile path/to/search_results.pep.xml \
  --outputDir ./results \
  --allele HLA-A*02:01 \
  --model Percolator \
  --testFDR 0.01
```

## 3. Inspect the Output

After a successful run, the output directory will contain:

```
results/
├── my_first_run.mokapot.psms.txt       # Rescored PSMs with q-values
├── my_first_run.mokapot.peptides.txt   # Peptide-level results
├── my_first_run.pin                     # PIN file with all features
├── models/                              # Saved rescoring model(s)
└── figures/
    ├── qvalues.png                      # Q-value curves
    ├── feature_importance.png           # Feature importance bar chart
    ├── feature_correlation.png          # Feature correlation heatmap
    └── target_decoy_histogram.png       # Target vs decoy distributions
```

## What's Next?

- See [Examples](../tutorial/examples.md) for complete Class I, Class II, and experiment-mode configurations
- Read [Pipeline Workflow](../tutorial/workflow.md) for a detailed explanation of each pipeline step
- Explore [Features](../tutorial/features/index.md) to understand what each feature computes
