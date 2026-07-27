# Configuration Reference

OptiMHC reads YAML configuration files. Values in the file override the defaults
defined by the application.

## Pipeline settings

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `experimentName` | String | `optimhc_experiment` | Output subdirectory name. |
| `inputType` | String | `pepxml` | Input format: `pepxml` or `pin`. |
| `inputFile` | String or list | Required | One or more PSM input files. |
| `decoyPrefix` | String | `DECOY_` | Protein prefix used to identify decoys in pepXML. |
| `outputDir` | String | `./results` | Base output directory. |
| `visualization` | Boolean | `true` | Generate result plots. |
| `saveModels` | Boolean | `true` | Save trained rescoring models. |
| `numProcesses` | Integer | `4` | Worker processes used by supported steps. |
| `showProgress` | Boolean | `true` | Display progress information. |
| `logLevel` | String | `INFO` | Logging level. |
| `keepIntermediate` | Boolean | `true` | Write supported intermediate prediction tables. |
| `toFlashLFQ` | Boolean | `true` | Write peptides passing `rescore.testFDR` in FlashLFQ format. |
| `retentionTimeColumn` | String | Auto-detected | PIN column containing retention time. |
| `allele` | List | `[]` | Alleles used by binding-affinity and PWM features. |
| `featureGenerator` | List | `[]` | Feature generators and their parameters. |
| `rescore` | Mapping | See below | Mokapot model and FDR settings. |

## Feature generators

Each generator is configured with a `name` and optional `params`:

```yaml
featureGenerator:
  - name: Basic
  - name: DeepLC
    params:
      calibrationCriteria: expect
      lowerIsBetter: true
      calibrationSize: 0.1
  - name: SpectralSimilarity
    params:
      mzmlDir: /path/to/mzml
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: LUMOS
      tolerance: 20
      numTopPeaks: 36
```

See the [feature reference](../tutorial/features/index.md) for generator-specific
requirements and parameters.

## Rescoring

```yaml
rescore:
  testFDR: 0.01
  trainFDR: 0.01
  model: Percolator
  numJobs: 4
  seed: 1
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `testFDR` | Float | `0.01` | FDR threshold used to evaluate and report results. |
| `trainFDR` | Float | `0.01` | FDR threshold used during model training. |
| `model` | String | `Percolator` | `Percolator`, `XGBoost`, or `RandomForest`. |
| `numJobs` | Integer | `1` | Parallel jobs used by supported models. |
| `seed` | Integer | `1` | Random seed passed to the model and Mokapot. |

See [Configuration Examples](../tutorial/examples.md) for complete Class I,
Class II, and experiment-mode files.
