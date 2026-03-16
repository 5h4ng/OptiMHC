import json
import logging

import click

from optimhc import __version__
from optimhc.core import Pipeline
from optimhc.core.config import Config

logger = logging.getLogger(__name__)

LOG_MAPPING = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def setup_logging(level: str = "INFO") -> None:
    if level not in LOG_MAPPING:
        raise ValueError(f"Invalid log level: {level}")
    logging.basicConfig(
        level=LOG_MAPPING[level],
        format="%(asctime)s [%(levelname)s] %(name)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )

    # mhctools attaches its own INFO-level handlers to its loggers
    # https://github.com/openvax/mhctools/blob/master/mhctools/logging.conf
    for name in [
        "mhctools",
        "mhctools.base_commandline_predictor",
        "mhctools.netmhc",
        "mhctools.netmhciipan",
        "mhctools.process_helpers",
        "mhctools.cleanup_context",
    ]:
        lg = logging.getLogger(name)
        lg.handlers.clear()
        lg.disabled = True
        lg.propagate = False
        lg.setLevel(logging.CRITICAL)


@click.group()
@click.version_option(version=__version__, prog_name="optimhc")
def cli():
    """
    OptiMHC - A optimized rescoring pipeline for immunopeptidomics data.
    """
    pass


@cli.command()
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to YAML configuration file",
)
@click.option(
    "--inputType",
    type=click.Choice(["pepxml", "pin"]),
    help="Type of input file",
)
@click.option(
    "--inputFile",
    type=click.Path(exists=True),
    multiple=True,
    help="Path(s) to input PSM file(s). Can be specified multiple times for multiple files.",
)
@click.option(
    "--decoyPrefix",
    type=str,
    help="Prefix used to identify decoy sequences",
)
@click.option(
    "--outputDir",
    type=click.Path(),
    help="Output directory",
)
@click.option(
    "--visualization/--no-visualization",
    is_flag=True,
    default=None,
    help="Enable/disable visualization",
)
@click.option(
    "--numProcesses",
    type=int,
    help="Number of parallel processes",
)
@click.option(
    "--allele",
    type=str,
    multiple=True,
    help="Allele(s) for which predictions will be computed",
)
@click.option(
    "--featureGenerator",
    type=str,
    multiple=True,
    help="Feature generator configuration in JSON format",
)
@click.option(
    "--logLevel",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"]),
    help="Logging level",
)
@click.option(
    "--testFDR",
    type=float,
    help="FDR threshold for testing",
)
@click.option(
    "--model",
    type=click.Choice(["Percolator", "XGBoost", "RandomForest"]),
    help="Model to use for rescoring",
)
def pipeline(
    config,
    inputtype,
    inputfile,
    decoyprefix,
    outputdir,
    visualization,
    numprocesses,
    allele,
    featuregenerator,
    loglevel,
    testfdr,
    model,
):
    """Run the optiMHC pipeline with the specified configuration."""
    pipeline_config = Config(config) if config else Config()

    if inputtype:
        pipeline_config["inputType"] = inputtype
    if inputfile:
        pipeline_config["inputFile"] = list(inputfile)
    if decoyprefix:
        pipeline_config["decoyPrefix"] = decoyprefix
    if outputdir:
        pipeline_config["outputDir"] = outputdir
    if visualization is not None:
        pipeline_config["visualization"] = visualization
    if numprocesses:
        pipeline_config["numProcesses"] = numprocesses
    if allele:
        pipeline_config["allele"] = list(allele)
    if loglevel:
        pipeline_config["logLevel"] = loglevel
    if featuregenerator:
        feature_generators = []
        for fg in featuregenerator:
            try:
                fg_config = json.loads(fg)
                feature_generators.append(fg_config)
            except json.JSONDecodeError as e:
                raise click.BadParameter(f"Invalid JSON format for feature generator: {e}")
        pipeline_config["featureGenerator"] = feature_generators
    if testfdr:
        pipeline_config["rescore"]["testFDR"] = testfdr
    if model:
        pipeline_config["rescore"]["model"] = model

    setup_logging(pipeline_config["logLevel"])
    pipeline_config.validate()
    Pipeline(pipeline_config).run()


@cli.command()
@click.option(
    "--config",
    type=click.Path(exists=True),
    required=True,
    help="Path to YAML configuration file",
)
def experiment(config):
    """Run multiple experiments with different feature combinations."""
    pipeline_config = Config(config)
    setup_logging(pipeline_config["logLevel"])

    Pipeline(pipeline_config).run_experiments()


if __name__ == "__main__":
    cli()
