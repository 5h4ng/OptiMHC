# OptiMHC

[<img src="img/optimhc_logo.svg" align="right" width="260" alt="OptiMHC logo">](https://github.com/5h4ng/OptiMHC)

**An optimized rescoring pipeline for immunopeptidomics data that significantly enhances peptide identification performance.**

<br clear="right"/>

## Quick Start

### Installation

<!-- ```bash
pip install optimhc
``` -->

```bash
git clone https://github.com/5h4ng/OptiMHC.git
cd OptiMHC
pip install -e .
```

## Usage

Create a YAML configuration:

```yaml
experimentName: example
inputType: pepxml
inputFile:
  - /path/to/results.pep.xml
outputDir: ./results

featureGenerator:
  - name: Basic

rescore:
  model: Percolator
```

```bash
optimhc pipeline --config /path/to/config.yaml
```

## Documentation

- [Quick Start](docs/getting-started/quickstart.md)
- [Configuration Reference](docs/getting-started/configuration.md)
- [Configuration Examples](docs/tutorial/examples.md)
- [Pipeline Workflow](docs/tutorial/workflow.md)
- [Feature Reference](docs/tutorial/features/index.md)
- [API Reference](https://optimhc.readthedocs.io/)

Release notes and migration information are available in the [changelog](CHANGELOG.md).
