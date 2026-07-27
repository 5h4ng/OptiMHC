# Examples

This page provides configuration templates that you can copy into local YAML files and adapt to your data. It explains each section of Class I, Class II, and experiment-mode configurations.

## MHC Class I Example

This configuration demonstrates a full Class I immunopeptidomics rescoring workflow with all available features.

```yaml
experimentName: classI_example
inputType: pepxml
inputFile:
  - /path/to/class_i_search_results.pep.xml
decoyPrefix: DECOY_
outputDir: ./results
visualization: True
numProcesses: 32
showProgress: True
keepIntermediate: True
toFlashLFQ: True
```

**General settings:**

- `experimentName` — the output subdirectory name under `outputDir`.
- `inputType` — the format of your search engine output (`pepxml` or `pin`).
- `inputFile` — one or more paths to PepXML files.
- `decoyPrefix` — the prefix used by the search engine to mark decoy protein accessions (default: `DECOY_`).
- `outputDir` — where results, models, and figures are written.
- `numProcesses` — number of parallel processes for feature generation.
- `keepIntermediate` — write supported intermediate results such as the best-prediction `BA.parquet` summary (default: `True`).
- `toFlashLFQ` — write peptides accepted at `rescore.testFDR` to a FlashLFQ-compatible table
  (default: `True`).

### Modification Handling

PIN and pepXML readers identify supported modifications using the bundled AlphaBase table.
DeepLC uses the parsed modification names and sites, while SpectralSimilarity receives ProForma
peptide strings.

### Allele Settings

```yaml
allele:
  - HLA-A*02:02
```

Specifies the MHC allele(s) for binding prediction tools. Use standard HLA nomenclature in config files.

The same allele list is used by NetMHCpan, NetMHCIIpan, MHCflurry, and PWM when those features are enabled. For Class I, use names such as `HLA-A*02:02` or `HLA-B*07:02`. For Class II, current examples use paired DP/DQ names such as `HLA-DPA1*02:01-DPB1*01:01`. PWM requires a matching local matrix under `optimhc/PWMs/`; see [PWM Score](features/pwm.md#allele-support-and-format).

### Features

```yaml
featureGenerator:
  - name: Basic
  - name: SpectralSimilarity
    params:
      mzmlDir: /path/to/mzml
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: LUMOS
      tolerance: 20
      numTopPeaks: 36
      url: 127.0.0.1:8500
      ssl: false
  - name: DeepLC
    params:
      calibrationCriteria: expect
      lowerIsBetter: True
      calibrationSize: 0.1
  - name: OverlappingPeptide
    params:
      minOverlapLength: 7
      minLength: 7
      maxLength: 20
      overlappingScore: expect
  - name: PWM
    params:
      class: I
  - name: MHCflurry
  - name: NetMHCpan
```

Each entry in `featureGenerator` specifies a feature `name` and optional `params`. Features without `params` use their defaults. See the [Features](features/index.md) section for detailed documentation of each feature.

### Rescoring Settings

```yaml
rescore:
  testFDR: 0.01
  model: Percolator
  numJobs: 4
```

- `testFDR` — the FDR threshold for the test set (default: 0.01).
- `model` — the rescoring model: `Percolator` (linear SVM), `XGBoost`, or `RandomForest`.
- `numJobs` — number of parallel jobs for cross-validation in XGBoost/RandomForest models.

---

## MHC Class II Example

This configuration mirrors the Class I example but is adapted for MHC Class II immunopeptidomics.

```yaml
experimentName: classII_example
inputType: pepxml
inputFile:
  - /path/to/class_ii_search_results.pep.xml
decoyPrefix: DECOY_
outputDir: ./results
```

### Key Differences from Class I

**Allele notation** uses the Class II alpha-beta chain format:

```yaml
allele:
  - HLA-DPA1*02:01-DPB1*01:01
```

**OverlappingPeptide** parameters are adjusted for the longer peptide lengths typical of Class II:

```yaml
- name: OverlappingPeptide
  params:
    minOverlapLength: 8
    minLength: 9
    maxLength: 50
```

**PWM** is set to Class II mode, which uses a sliding 9-mer core window with N- and C-flank scoring:

```yaml
- name: PWM
  params:
    class: II
```

**NetMHCIIpan** replaces MHCflurry and NetMHCpan (which are Class I only):

```yaml
- name: NetMHCIIpan
```

---

## Experiment Mode Example

Experiment mode runs multiple rescoring experiments with different feature subsets on the same input data, allowing you to compare the contribution of individual features.

The general settings and features are defined once at the top. The `experiments` section then defines each experiment:

```yaml
experiments:
  - name: "Baseline"
    source: ["Original"]
    model: "Percolator"
  - name: "Complete"
    source: ["Original", "OverlappingPeptide", "ContigFeatures", "PWM", "Basic"]
    model: "Percolator"
  - name: "Shuffle"
    source: ["Original", "Basic", "OverlappingPeptide", "ContigFeatures", "PWM"]
    model: "Percolator"
```

Each experiment specifies:

- `name` — a label for the experiment, used for output subdirectories.
- `source` — complete feature categories. `Original` contains reader-provided search features;
  each configured generator contributes a category with the same name.
- `features` — optional advanced alternative for selecting exact DataFrame columns. Do not use
  `source` and `features` in the same experiment.
- `model` — the rescoring model to use for this experiment.

Save the completed configuration as a local YAML file, then run experiment mode with:

```bash
optimhc experiment --config experiment_example.yaml
```

!!! tip
    Prefer `source` for stable experiment YAML. Generator implementations may add or rename
    columns without requiring the YAML to enumerate every output feature.
