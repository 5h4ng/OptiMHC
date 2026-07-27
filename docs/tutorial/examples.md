# Examples

These templates enable the full feature set for HLA-I or HLA-II. Copy one into
a YAML file, replace the input paths and alleles, and check the external tool
requirements before running it.

Default values are omitted unless they make the template easier to understand.
See the [Configuration Reference](../getting-started/configuration.md) for all
settings and defaults.

## HLA-I full model

```yaml
experimentName: hla_i_full_model
inputType: pepxml
inputFile:
  - /data/run_01.pep.xml
  - /data/run_02.pep.xml
outputDir: ./results
numProcesses: 36

allele:
  - HLA-A*02:01
  - HLA-B*07:02
  - HLA-C*07:02

featureGenerator:
  - name: Basic
  - name: SpectralSimilarity
    params:
      mzmlDir: /data/mzml
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: QE
      tolerance: 20
      numTopPeaks: 36
      url: koina.wilhelmlab.org:443
  - name: DeepLC
    params:
      calibrationCriteria: expect
      lowerIsBetter: true
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

rescore:
  model: Percolator
```

This writes results to `./results/hla_i_full_model/`. `MHCflurry` is installed
with OptiMHC. `NetMHCpan` must be installed separately and available on `PATH`.

## HLA-II full model

```yaml
experimentName: hla_ii_full_model
inputType: pepxml
inputFile:
  - /data/run_01.pep.xml
  - /data/run_02.pep.xml
outputDir: ./results
numProcesses: 36

allele:
  - DRB1*15:01
  - HLA-DPA1*02:01-DPB1*01:01
  - HLA-DQA1*05:01-DQB1*02:01

featureGenerator:
  - name: Basic
  - name: SpectralSimilarity
    params:
      mzmlDir: /data/mzml
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: QE
      tolerance: 20
      numTopPeaks: 36
      url: koina.wilhelmlab.org:443
  - name: DeepLC
    params:
      calibrationCriteria: expect
      lowerIsBetter: true
  - name: OverlappingPeptide
    params:
      minOverlapLength: 8
      minLength: 9
      maxLength: 50
      overlappingScore: expect
  - name: PWM
    params:
      class: II
  - name: NetMHCIIpan

rescore:
  model: Percolator
```

This writes results to `./results/hla_ii_full_model/`. `NetMHCIIpan` must be
installed separately and available on `PATH`.

## Before running a full model

- Replace both `inputFile` paths and `mzmlDir`.
- Confirm that each PSM run maps to `<mzmlDir>/<run>.mzML`.
- Confirm that `expect` exists and lower values are better.
- Use alleles supported by every configured predictor and by the bundled PWM files.
- Add `params.executablePath` if NetMHCpan or NetMHCIIpan is not on `PATH`.

`SpectralSimilarity` sends prediction inputs to the configured Koina server.
The templates use the public endpoint with SSL. Use a self-hosted endpoint when
the prediction inputs must remain local.

## Experiment mode

Experiment mode compares feature subsets on the same input data. General
settings and generators remain at the top level. The `experiments` list selects
the feature groups and model for each run.

```yaml
experiments:
  - name: Baseline
    source: [Original]
    model: Percolator
  - name: Complete
    source: [Original, Basic, OverlappingPeptide, ContigFeatures, PWM]
    model: Percolator
```

```bash
optimhc experiment --config experiment_example.yaml
```

!!! tip
    Prefer `source` when selecting complete feature groups. Use `features` only
    when an experiment needs exact DataFrame columns. Do not set both.
