import importlib.util
import json
import logging
import os
import sys

import click

from optimhc.core import Pipeline
from optimhc.core.config import Config

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[logging.StreamHandler()],
)

logger = logging.getLogger(__name__)


@click.group()
def cli():
    """
    optiMHC - A high-performance rescoring pipeline for immunopeptidomics data.
    """
    pass


def parse_cli_config(**kwargs):
    # Remove None values and build a config dict
    return {k: v for k, v in kwargs.items() if v is not None and v != ()}


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
    # Load configuration
    if config:
        pipeline_config = Config(config)
    else:
        pipeline_config = Config()

    # Override with command-line parameters
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
        pipeline_config["numProcess"] = numprocesses
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

    # Run pipeline
    pipeline_config.validate()
    pipeline = Pipeline(pipeline_config)
    pipeline.run()


@cli.command()
@click.option(
    "--config",
    type=click.Path(exists=True),
    required=True,
    help="Path to YAML configuration file",
)
def experiment(config):
    """Run multiple experiments with different feature combinations."""
    # Load configuration
    pipeline_config = Config(config)

    # Run experiments
    pipeline = Pipeline(pipeline_config)
    pipeline.run_experiments()


@cli.command()
def gui():
    """Launch the optiMHC GUI."""
    if importlib.util.find_spec("streamlit") is None:
        print("Error: Streamlit is not installed. Install GUI dependencies with:")
        print("pip install optimhc[gui]")
        return

    import subprocess

    # Get the path to the GUI app
    gui_path = os.path.join(os.path.dirname(__file__), "gui", "app.py")

    if not os.path.exists(gui_path):
        print(f"Error: GUI application not found at {gui_path}")
        return

    # Create a temporary launcher script that uses the correct imports
    import tempfile

    launcher_content = """
import os
import sys
import streamlit

# Add the root directory to the path
sys.path.insert(0, '{}')

# Import the app module properly
from optimhc.gui.app import main

if __name__ == "__main__":
    main()
    """.format(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

    fd, temp_path = tempfile.mkstemp(suffix=".py")
    with os.fdopen(fd, "w") as f:
        f.write(launcher_content)

    # Launch Streamlit with the temporary script
    print("Starting optiMHC GUI...")
    try:
        subprocess.run([sys.executable, "-m", "streamlit", "run", temp_path])
    finally:
        # Clean up the temporary file
        try:
            os.unlink(temp_path)
        except OSError:
            pass


if __name__ == "__main__":
    cli()
