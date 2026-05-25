# Examples

OptiMHC ships with example configuration files under the `examples/` directory. This page walks through each one, explaining every section.

## MHC Class I Example

**File:** `examples/classI_example.yaml`

This configuration demonstrates a full Class I immunopeptidomics rescoring workflow with all available features.

```yaml
experimentName: classI_example
inputType: pepxml
inputFile:
  - ../data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml
decoyPrefix: DECOY_
outputDir: ./examples/results
visualization: True
removePreNxtAA: False
numProcesses: 32
showProgress: True
```

**General settings:**

- `experimentName` — a label used for output file naming.
- `inputType` — the format of your search engine output (`pepxml` or `pin`).
- `inputFile` — one or more paths to PepXML files.
- `decoyPrefix` — the prefix used by the search engine to mark decoy protein accessions (default: `DECOY_`).
- `outputDir` — where results, models, and figures are written.
- `numProcesses` — number of parallel processes for feature generation.

### Modification Mapping

```yaml
modificationMap:
  "147.035385": "UNIMOD:35"  # Oxidation (M)
  "160.030649": "UNIMOD:4"   # Carbamidomethyl (C)
```

Maps the **full modified residue mass** (residue + modification, as recorded in PepXML) to a UNIMOD identifier. This is required by generators that need standardized modification notation (e.g., SpectralSimilarity, DeepLC).

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
      mzmlDir: ../data
      spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+
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

**File:** `examples/classII_example.yaml`

This configuration mirrors the Class I example but is adapted for MHC Class II immunopeptidomics.

```yaml
experimentName: classII_example
inputType: pepxml
inputFile:
  - ../data/AG20201214_FAIMS_DPB0101_DPA0201_93e6_1hr.pep.xml
decoyPrefix: DECOY_
outputDir: ./examples/results
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

**File:** `examples/experiment_example.yaml`

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
- `source` — a list of feature sources to include. These correspond to the source names registered by each feature (e.g., `"Original"` from the parser, `"Basic"` from Basic, etc.).
- `model` — the rescoring model to use for this experiment.

Run experiment mode with:

```bash
optimhc experiment --config examples/experiment_example.yaml
```

!!! tip
    Experiment mode is useful for ablation studies — start with `"Original"` as a baseline, then incrementally add feature sources to measure their impact on peptide identification.

### Available Feature Sources

| Source Name | Feature |
|---|---|
| `Original` | PepXML / PIN parser (search engine features) |
| `Basic` | Basic peptide features |
| `SpectralSimilarity` | Spectral similarity |
| `DeepLC` | Retention time deviation |
| `OverlappingPeptide` | Overlapping peptide score |
| `ContigFeatures` | Contig-level features from overlapping peptide analysis |
| `PWM` | Position weight matrix score |
| `MHCflurry` | MHCflurry binding, processing, and presentation prediction |
| `NetMHCpan` | NetMHCpan class I binding affinity prediction |
| `NetMHCIIpan` | NetMHCIIpan class II binding affinity prediction |
