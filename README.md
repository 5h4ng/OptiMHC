# optiMHC

**optiMHC** is a rescoring pipeline for immunopeptidomics data. It enhances peptide identification by integrating multiple feature generators and machine learning-based rescoring. The package is modular, supports extensible feature generation, and provides both command-line and programmatic interfaces.

## Quick Start

> TODO

## Usage

### Using a YAML Configuration File (Recommended)

Using a YAML configuration file is recommended because it provides a more flexible and user-friendly way to configure the pipeline.

```bash
optimhc pipeline --config /path/to/config.yaml
```

**Note:** The default configuration is stored in `optimhc/core/config.py`. Your custom configuration will override the default values.

#### Configuration Parameters

The pipeline can be configured by using a YAML file. This file defines the input settings, the list of feature generators, rescore parameters, and (optionally) experiment configurations. Below you will find a table summarizing the main configuration parameters along with examples and descriptions.

| Parameter           | Type               | Example                                                          | Description                                                             |
|---------------------|--------------------|------------------------------------------------------------------|-------------------------------------------------------------------------|
| `experimentName`    | String             | `classI_example`                                                 | Name of the experiment and output subdirectory name.                  |
| `inputType`         | String             | `pepxml`                                                         | Type of input file. Supported values: `pepxml`, `pin`.                  |
| `inputFile`         | String or List     | `./data/xxx.pep.xml`        | Path(s) to the input PSM file(s).                                       |
| `decoyPrefix`       | String             | `DECOY_`                                                         | Prefix used to identify decoy sequences.                              |
| `outputDir`         | String             | `./results`                                                      | Base directory where output files, logs and figures are stored.         |
| `visualization`     | Boolean            | `True`                                                           | Enable or disable generation of visualization plots.                  |
| `removePreNxtAA`    | Boolean            | `False`                                                          | Remove pre/post neighboring amino acids in sequence processing.         |
| `numProcesses`      | Integer            | `32`                                                             | Number of parallel processes to use.                               |
| `showProgress`      | Boolean            | `True`                                                           | Show progress information during execution.     
| `logLevel`          | String             | `INFO`                                                           | Logging level (DEBUG, INFO, WARNING, ERROR). Default is "INFO".
| `modificationMap`   | Dictionary         | `{ '147.035385': 'UNIMOD:35' }`                                    | Maps modification masses to their 'UNIMOD' identifiers. See https://www.unimod.org/ for details           |
| `allele`            | List               | `[HLA-A*02:02]`                                                  | List of alleles for which predictions will be computed.                 |
| `featureGenerator`  | List of Dictionaries | See table below                                                  | List of feature generator configurations (each with a `name` and optional `params`). |
| `rescore`           | Dictionary         | See table below                                                 | Rescore settings including FDR threshold, model and number of jobs.     |

#### Feature Generator Configurations

Each feature generator is specified with its `name` and an optional `params` subsection. Some common generators include:

| Generator Name      | Example Parameters                                                                                                           | Description                                                                              |
|---------------------|------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------|
| `Basic`             | N/A                                                                                                                          | Generates basic sequence features.                                                       |
| `SpectraSimilarity` | `mzmlDir: ./data`<br>`spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+`<br>`model: AlphaPeptDeep_ms2_generic`<br>`collisionEnergy: 28`<br>`instrument: LUMOS`<br>`tolerance: 20`<br>`numTopPeaks: 36`<br>`url: koina.wilhelmlab.org:443` | Computes features based on the similarity between experimental spectra and predicted spectra. See more options on https://koina.proteomicsdb.org/ |
| `DeepLC`            | `calibrationCriteria: expect`<br>`lowerIsBetter: True`<br>`calibrationSize: 0.1`                                             | Creates retention time predictions by calibrating using DeepLC.                           |
| `OverlappingPeptide`     | `minOverlapLength: 7`<br>`minLength: 7`<br>`maxLength: 20`<br>`overlappingScore: expect`                                           | Generates overlapping peptide features for grouping similar peptides.                         |
| `PWM`               | `class: I`                                                                                                                   | Generates position weight matrix features for MHC class I and class II peptides.                        |
| `MHCflurry`         | N/A                                                                                                                          | Predicts class I binding affinities using the MHCflurry model.                                     |
| `NetMHCpan`         | N/A                                                                                                                          | Predicts class I peptide-MHC binding affinity using NetMHCpan.                                     |
| `NetMHCIIpan` | N/A                                                                                                                          | Predicts class II peptide-MHC binding affinity using NetMHCIIpan.                                   |

#### Rescore Settings

Rescore parameters control how the rescoring step is executed and include:

| Parameter  | Type    | Example     | Description                                                                |
|------------|---------|-------------|----------------------------------------------------------------------------|
| `testFDR`  | Float   | `0.01`      | The false-discovery rate threshold at which to evaluate the learned models.       |
| `model`    | String  | `Percolator`| Model to use for rescoring (valid options include `Percolator`, `XGBoost`, or `RandomForest`). |
| `numJobs`  | Integer | `4`         |The number of parallel jobs to run. This value is passed to Scikit-learn's n_jobs parameter to control parallelism for model training or scoring. Set to -1 to use all available CPU cores.                   |

#### Example YAML Configuration

Below is an example YAML configuration for class I based on the latest pipeline version:

```yaml
experimentName: classI_example
inputType: pepxml
inputFile:
  - ./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml
decoyPrefix: DECOY_
outputDir: ./results
visualization: True
removePreNxtAA: False
numProcesses: 32
showProgress: True
modificationMap:
  '147.035385': 'UNIMOD:35'

# Allele settings
allele:
  - HLA-A*02:02

# Feature generator configurations
featureGenerator:
  - name: Basic
  - name: SpectraSimilarity
    params:
      mzmlDir: ./data
      spectrumIdPattern: (.+?)\.\d+\.\d+\.\d+
      model: AlphaPeptDeep_ms2_generic
      collisionEnergy: 28
      instrument: LUMOS
      tolerance: 20  
      numTopPeaks: 36
      url: koina.wilhelmlab.org:443
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

# Rescore settings
rescore:
  testFDR: 0.01
  model: Percolator
  numJobs: 4
```

### Using Direct Command-Line Parameters (Optional)

While we recommend using the YAML configuration file, you can also use command-line parameters to configure the pipeline:

```bash
optimhc pipeline \
  --inputType pepxml \
  --inputFile ./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml \
  --decoyPrefix DECOY_ \
  --outputDir ./results \
  --visualization \
  --numProcesses 32 \
  --allele HLA-A*02:02 \
  --logLevel INFO \
  --featureGenerator '{"name": "Basic"}' \
  --testFDR 0.01 \
  --model Percolator
```

**Note:** If you use both YAML configuration file and command-line parameters, command-line parameters will override the corresponding values in the YAML configuration file. 

#### Feature Generator Command-line Parameters 

The `--featureGenerator` option accepts JSON formatted strings that define the feature generator configuration. You can specify multiple feature generators by using the option multiple times. 

But be careful that if you use `--featureGenerator` in command-line, all your feature generator configurations in YAML file (`--config`) will be ignored.

Thus, **rather than using both methods simultaneously, use either command-line arguments or YAML for feature generator configuration.**

Here are some examples:

<details>
<summary>Basic feature generator (no parameters)</summary>

```bash
--featureGenerator '{"name": "Basic"}'
```

</details>

<details>
<summary>SpectraSimilarity with parameters</summary>

```bash
--featureGenerator '{
  "name": "SpectraSimilarity",
  "params": {
    "mzmlDir": "./data",
    "spectrumIdPattern": "(.+?)\.\d+\.\d+\.\d+",
    "model": "AlphaPeptDeep_ms2_generic",
    "collisionEnergy": 28,
    "instrument": "LUMOS",
    "tolerance": 20,
    "numTopPeaks": 36,
    "url": "koina.wilhelmlab.org:443"
  }
}'
```

</details>

<details>
<summary>Multiple feature generators</summary>

```bash
--featureGenerator '{"name": "Basic"}' \
--featureGenerator '{
  "name": "SpectraSimilarity",
  "params": {
    "mzmlDir": "./data",
    "model": "AlphaPeptDeep_ms2_generic"
  }
}' \
--featureGenerator '{
  "name": "DeepLC",
  "params": {
    "calibrationCriteria": "expect",
    "lowerIsBetter": true,
    "calibrationSize": 0.1
  }
}'
```

</details>

<details>
<summary>Some tips for JSON format</summary>

- Use single quotes (`'`) to wrap the entire JSON string
- All JSON strings must be valid JSON format (e.g., use `true` instead of `True`, `false` instead of `False`)
- For complex parameters, you can use a single line with proper escaping:
```bash
--featureGenerator '{"name":"SpectraSimilarity","params":{"mzmlDir":"./data","model":"AlphaPeptDeep_ms2_generic"}}'
```

</details>

### Full CLI Help

```bash
optimhc --help
optimhc pipeline --help
optimhc experiment --help
```
