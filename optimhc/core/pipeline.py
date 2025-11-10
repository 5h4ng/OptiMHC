"""
pipeline.py

Implements the main optiMHC pipeline for immunopeptidomics rescoring, including input parsing,
feature generation, rescoring, result saving, and visualization. Supports both single-run and
experiment modes.
"""
import os
import logging
import gc
import pandas as pd
from multiprocessing import Process
from mokapot.model import PercolatorModel

from optimhc.parser import read_pin, read_pepxml
from optimhc.rescore import mokapot
from optimhc.rescore.model import XGBoostPercolatorModel, RandomForestPercolatorModel
from optimhc.visualization import (
    plot_feature_importance,
    visualize_target_decoy_features,
    visualize_feature_correlation,
    plot_qvalues,
)
from optimhc.core.config import load_config, Config
from optimhc.core.logging_helper import setup_loggers
from optimhc.core.feature_generation import generate_features

logger = logging.getLogger(__name__)

# TODO: fix the type hinting for mokapot results
class Pipeline:
    """
    Main pipeline class for optiMHC, encapsulating the full data processing workflow.

    This class orchestrates input parsing, feature generation, rescoring, result saving, and visualization.
    It supports both single-run and experiment modes (multiple feature/model combinations).

    Parameters
    ----------
    config : str, dict, or Config
        Path to YAML config, dict, or Config object.

    Examples
    --------
    >>> from optimhc.core import Pipeline
    >>> pipeline = Pipeline(config)
    >>> pipeline.run()
    """
    def __init__(self, config):
        """
        Initialize the pipeline with a configuration file, dict, or Config object.

        Parameters
        ----------
        config : str, dict, or Config
            Path to YAML config, dict, or Config object.
        """
        logger.debug(f"config: {config}")
        if isinstance(config, Config):
            self.config = config
        else:
            self.config = Config(config)
        self.config.validate()
        self.experiment = self.config.get("experimentName", "optimhc_experiment")
        self.output_dir = os.path.join(self.config["outputDir"], self.experiment)
        os.makedirs(self.output_dir, exist_ok=True)
        setup_loggers(os.path.join(self.output_dir, "log"), self.config.get("logLevel", "INFO"))

        self.visualization_enabled = self.config.get("visualization", True)
        self.save_models = self.config.get("saveModels", True)
        self.to_flashlfq = self.config.get("toFlashLFQ", True)
        self.test_fdr = self.config.get("rescore", {}).get("testFDR", 0.01)
        self.model_type = self.config.get("rescore", {}).get("model", "Percolator")
        self.n_jobs = self.config.get("rescore", {}).get("numJobs", 1)

    def read_input(self):
        """
        Read input PSMs based on configuration.

        Returns
        -------
        PsmContainer
            Object containing loaded PSMs.

        Raises
        ------
        ValueError
            If input type is unsupported.
        Exception
            If file reading fails.
        """
        input_type = self.config["inputType"]
        input_files = self.config["inputFile"]
        if not isinstance(input_files, list):
            input_files = [input_files]

        try:
            if input_type == "pepxml":
                psms = read_pepxml(input_files, decoy_prefix=self.config["decoyPrefix"])
            elif input_type == "pin":
                psms = read_pin(input_files, retention_time_column=self.config.get("retentionTimeColumn"))
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            return psms
        except Exception as e:
            logger.error(f"Failed to read input files: {e}")
            raise

    def _generate_features(self, psms):
        """
        Generate features for PSMs using the configured feature generators.

        Parameters
        ----------
        psms : PsmContainer
            PSM container object.

        Returns
        -------
        PsmContainer
            PSM container with generated features.
        """
        generate_features(psms, self.config)
        return psms

    def rescore(self, psms, model_type=None, n_jobs=None, test_fdr=None, rescoring_features=None):
        """
        Perform rescoring on the PSMs using the specified or configured model.

        Parameters
        ----------
        psms : PsmContainer
            PSM container object.
        model_type : str, optional
            Model type ('XGBoost', 'RandomForest', 'Percolator').
        n_jobs : int, optional
            Number of parallel jobs.
        test_fdr : float, optional
            FDR threshold.
        rescoring_features : list, optional
            List of features to use for rescoring.

        Returns
        -------
        results : mokapot.Results
            Rescoring results.
        models : list
            Trained models.

        Notes
        -----
        Rescoring logic is adapted from mokapot (https://mokapot.readthedocs.io/)
        """
        test_fdr = test_fdr if test_fdr is not None else self.test_fdr
        model_type = model_type if model_type is not None else self.model_type
        n_jobs = n_jobs if n_jobs is not None else self.n_jobs

        if model_type == "XGBoost":
            model = XGBoostPercolatorModel(n_jobs=n_jobs)
        elif model_type == "RandomForest":
            model = RandomForestPercolatorModel(n_jobs=n_jobs)
        elif model_type == "Percolator":
            model = PercolatorModel(n_jobs=n_jobs)
        else:
            model = PercolatorModel(n_jobs=n_jobs)

        kwargs = {}
        if rescoring_features is not None:
            kwargs["rescoring_features"] = rescoring_features

        results, models = mokapot.rescore(
            psms, model=model, test_fdr=test_fdr, **kwargs
        )
        return results, models

    def save_results(self, psms, results, models, output_dir=None, file_root="optimhc"):
        """
        Save rescoring results, PSM data, and trained models to disk.

        Parameters
        ----------
        psms : PsmContainer
            PSM container object.
        results : mokapot.Results
            Rescoring results.
        models : list
            Trained models.
        output_dir : str, optional
            Output directory.
        file_root : str, optional
            Root name for output files.
        """
        output_dir = output_dir if output_dir is not None else self.output_dir

        results.to_txt(dest_dir=output_dir, file_root=file_root, decoys=True)
        psms.write_pin(os.path.join(output_dir, f"{file_root}.pin"))

        if self.save_models:
            model_dir = os.path.join(output_dir, "models")
            os.makedirs(model_dir, exist_ok=True)
            logger.info(f"Saving models to {model_dir}")
            for i, model in enumerate(models):
                model.save(os.path.join(model_dir, f"{file_root}.model{i}"))
        
        if self.to_flashlfq:
            mokapot.to_flashLFQ(results, output_dir, file_name=f"{file_root}.FlashLFQ.txt")

    def visualize_results(self, psms, results, models, output_dir=None, sources=None):
        """
        Generate and save visualizations for the analysis results.

        Parameters
        ----------
        psms : PsmContainer
            PSM container object.
        results : mokapot.Results
            Rescoring results.
        models : list
            Trained models.
        output_dir : str, optional
            Output directory.
        sources : list, optional
            Feature sources to include in visualizations.
        """
        if not self.visualization_enabled:
            logger.info("Visualization is disabled. Skipping...")
            return

        output_dir = output_dir if output_dir is not None else self.output_dir
        fig_dir = os.path.join(output_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        plot_qvalues(
            results,
            save_path=os.path.join(fig_dir, "qvalues.png"),
            threshold=0.05,
        )

        if sources:
            rescoring_features = {
                k: v for k, v in psms.rescoring_features.items() if k in sources
            }
        else:
            rescoring_features = psms.rescoring_features

        plot_feature_importance(
            models,
            rescoring_features,
            save_path=os.path.join(fig_dir, "feature_importance.png"),
        )
        visualize_feature_correlation(
            psms,
            save_path=os.path.join(fig_dir, "feature_correlation.png"),
        )
        visualize_target_decoy_features(
            psms,
            num_cols=4,
            save_path=os.path.join(fig_dir, "target_decoy_histogram.png"),
        )

    def _run_single_experiment(self, psms, exp_config, exp_name, exp_dir):
        """
        Run a single experiment with the specified configuration.

        Parameters
        ----------
        psms : PsmContainer
            PSM container object.
        exp_config : dict
            Experiment-specific configuration.
        exp_name : str
            Name of the experiment.
        exp_dir : str
            Output directory for the experiment.

        Returns
        -------
        bool
            True if experiment succeeded, False otherwise.
        """
        try:
            os.makedirs(exp_dir, exist_ok=True)

            source = exp_config.get("source", None)
            model_type = exp_config.get("model", self.model_type)
            n_jobs = exp_config.get("numJobs", self.n_jobs)

            logger.info(f"Running experiment '{exp_name}' with sources: {source}")

            # Generate list of features based on the provided sources
            features = []
            if source:
                for s in source:
                    if s in psms.rescoring_features:
                        features.extend(psms.rescoring_features.get(s, []))
                    else:
                        logger.error(f"Source '{s}' not found in PSM features")
            logger.info(f"Features used in experiment '{exp_name}': {features}")

            results, models = self.rescore(
                psms,
                model_type=model_type,
                n_jobs=n_jobs,
                test_fdr=self.test_fdr,
                rescoring_features=features,
            )

            self.save_results(
                psms, results, models, output_dir=exp_dir, file_root=exp_name
            )

            fig_dir = os.path.join(exp_dir, "figures")

            plot_qvalues(
                results,
                save_path=os.path.join(fig_dir, "qvalues.png"),
                threshold=0.05,
            )

            plot_feature_importance(
                models,
                rescoring_features={
                    k: v for k, v in psms.rescoring_features.items() if k in source
                },
                save_path=os.path.join(fig_dir, "feature_importance.png"),
            )

            return True

        except Exception as e:
            logger.error(f"Experiment '{exp_name}' failed: {e}")
            return False

        finally:
            # Explicit resource release to free up memory after each experiment
            try:
                del results
                del models
            except Exception:
                pass
            gc.collect()

    def run(self):
        """
        Run the complete optiMHC pipeline (single run mode).

        This method executes the full workflow: input parsing, feature generation, rescoring, saving, and visualization.

        Returns
        -------
        psms : PsmContainer
            PSM container object.
        results : mokapot.Results
            Rescoring results.
        models : list
            Trained models.
        """
        logger.info("Starting analysis pipeline")

        psms = self.read_input()
        psms = self._generate_features(psms)
        results, models = self.rescore(psms)
        self.save_results(psms, results, models)
        self.visualize_results(psms, results, models)

        logger.info(f"Analysis pipeline completed, results saved to {self.output_dir}")
        return psms, results, models

    def run_experiments(self):
        """
        Run experiments with different feature/model combinations using multiprocessing.

        Each experiment is executed in its own process for complete resource isolation.
        The experiment configurations must be provided in the config under the 'experiments' key.

        Returns
        -------
        None
        """
        logger.info("Starting experiment mode with multiple feature combinations")

        psms = self.read_input()
        psms = self._generate_features(psms)

        # Save the generated pin file for reference
        pin_path = os.path.join(self.output_dir, f"optimhc.{self.experiment}.pin")
        psms.write_pin(pin_path)
        fig_summary_dir = os.path.join(self.output_dir, "figures")
        os.makedirs(fig_summary_dir, exist_ok=True)
        visualize_feature_correlation(
            psms,
            save_path=os.path.join(fig_summary_dir, "feature_correlation.png"),
        )
        # visualize_target_decoy_features(
        #     psms,
        #     num_cols=4,
        #     save_path=os.path.join(fig_summary_dir, 'target_decoy_histogram.png'),
        # )

        experiment_configs = self.config.get("experiments", [])
        processes = []
        for i, exp_config in enumerate(experiment_configs):
            exp_name = exp_config.get("name", f"Experiment_{i + 1}")
            exp_dir = os.path.join(self.output_dir, exp_name)

            logger.info(f"Starting experiment '{exp_name}' in a separate process")
            p = Process(
                target=self._run_single_experiment,
                args=(psms, exp_config, exp_name, exp_dir),
            )
            p.start()
            processes.append(p)

        # Wait for all experiment processes to finish
        for p in processes:
            p.join()

        logger.info("All experiments completed")
