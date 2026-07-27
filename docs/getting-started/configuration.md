# Configuration Reference

OptiMHC combines three configuration sources in this order:

1. Command-line options
2. Values in the YAML file
3. Built-in defaults

You only need to write values that you want to change. Nested mappings such as
`rescore` are merged with their defaults. Lists such as `inputFile`, `allele`,
and `featureGenerator` replace the default list.

Run a YAML configuration with:

```bash
optimhc pipeline --config config.yaml
```

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
| `featureGenerator` | List | Not set | Feature generators and their parameters. |
| `rescore` | Mapping | See below | Mokapot model and FDR settings. |

The full templates set `numProcesses: 36` as a practical starting point for a
large workstation. The built-in default remains `4`. Reduce this value when
memory or CPU resources are limited.

## Input and output

`inputFile` accepts one path or a list of paths. All files in one run must use
the format selected by `inputType`:

```yaml
inputType: pepxml
inputFile:
  - /data/run_01.pep.xml
  - /data/run_02.pep.xml
```

`outputDir` is the base directory. `experimentName` creates one directory
inside it:

```yaml
experimentName: hla_i_full_model
outputDir: ./results
```

This example writes results to `./results/hla_i_full_model/`.

For `SpectralSimilarity`, each PSM run must have a matching mzML file at
`<mzmlDir>/<run>.mzML`.

## Alleles

Define all alleles once at the top level. The configured PWM and
binding-affinity generators reuse this list.

```yaml
allele:
  - HLA-A*02:01
  - HLA-B*07:02
  - HLA-C*07:02
```

Use standard HLA names. PWM also requires a matching bundled matrix. External
predictors validate alleles against their own supported lists.

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
      instrument: QE
      tolerance: 20
      numTopPeaks: 36
      url: koina.wilhelmlab.org:443
```

The public Koina server uses SSL by default, so `ssl: true` may be omitted.
Using `SpectralSimilarity` sends peptide prediction inputs to the configured
server. Set `url` to a self-hosted Koina endpoint if the data must remain local.

`NetMHCpan` and `NetMHCIIpan` use executables on `PATH`. If necessary, set
`params.executablePath` to the installed executable.

See the [feature reference](../tutorial/features/index.md) for generator-specific
requirements and parameters.

## Rescoring

```yaml
rescore:
  model: Percolator
```

| Parameter | Type | Default | Description |
| --- | --- | --- | --- |
| `testFDR` | Float | `0.01` | FDR threshold used to evaluate and report results. |
| `trainFDR` | Float | `0.01` | FDR threshold used during model training. |
| `model` | String | `Percolator` | `Percolator`, `XGBoost`, or `RandomForest`. |
| `numJobs` | Integer | `1` | Parallel jobs used by supported models. |
| `seed` | Integer | `1` | Random seed passed to the model and Mokapot. |

The example only states the model because all other values use their defaults.

See [Configuration Examples](../tutorial/examples.md) for complete HLA-I,
HLA-II, and experiment-mode files.
