import gc
import logging
import os
import re
from multiprocessing import Process

import pandas as pd

from optimhc.core.config import Config
from optimhc.core.feature_generation import generate_features, select_feature_groups
from optimhc.parser import read_pepxml, read_pin
from optimhc.rescore.factory import rescore_model_factory

logger = logging.getLogger(__name__)
CHARGE_FEATURE = re.compile(r"^charge_?(\d+)(?:_or_more)?$", re.IGNORECASE)


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

        self.visualization_enabled = self.config.get("visualization", True)
        self.save_models = self.config.get("saveModels", True)
        self.to_flashlfq = self.config.get("toFlashLFQ", True)
        self.test_fdr = self.config.get("rescore", {}).get("testFDR", 0.01)
        self.train_fdr = self.config.get("rescore", {}).get("trainFDR", 0.01)
        self.model_type = self.config.get("rescore", {}).get("model", "Percolator")
        self.n_jobs = self.config.get("rescore", {}).get("numJobs", 1)
        self.rng = self.config.get("rescore", {}).get("seed", 1)
        self.feature_groups = {}

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
                containers = [
                    read_pepxml(path, decoy_prefix=self.config["decoyPrefix"])
                    for path in input_files
                ]
            elif input_type == "pin":
                containers = [
                    read_pin(
                        path,
                        retention_time_column=self.config.get("retentionTimeColumn"),
                    )
                    for path in input_files
                ]
            else:
                raise ValueError(f"Unsupported input type: {input_type}")
            non_charge_features = tuple(
                feature
                for feature in containers[0].feature_columns
                if CHARGE_FEATURE.fullmatch(feature) is None
            )
            if any(
                tuple(
                    feature
                    for feature in container.feature_columns
                    if CHARGE_FEATURE.fullmatch(feature) is None
                )
                != non_charge_features
                for container in containers[1:]
            ):
                raise ValueError("All input files must declare the same rescoring features.")
            charge_features = sorted(
                {
                    feature
                    for container in containers
                    for feature in container.feature_columns
                    if CHARGE_FEATURE.fullmatch(feature) is not None
                },
                key=lambda feature: (
                    int(CHARGE_FEATURE.fullmatch(feature).group(1)),
                    feature.lower(),
                ),
            )
            feature_columns = (*non_charge_features, *charge_features)
            frames = []
            for container in containers:
                frame = container.df.copy()
                for feature in charge_features:
                    if feature not in frame:
                        frame[feature] = 0.0
                frames.append(frame)
            combined = pd.concat(frames, ignore_index=True)
            combined["psm_id"] = range(len(combined))
            return type(containers[0])(combined, feature_columns=feature_columns)
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
        generated = generate_features(psms, self.config)
        self.feature_groups = generated.feature_groups
        if self.config.get("keepIntermediate", True) and generated.raw_predictions:
            from optimhc.core.feature_generation import _build_ba_parquet

            intermediate_results_dir = os.path.join(self.output_dir, "intermediate_results")
            os.makedirs(intermediate_results_dir, exist_ok=True)
            _build_ba_parquet(
                generated.raw_predictions,
                os.path.join(intermediate_results_dir, "BA.parquet"),
            )
        return psms

    @staticmethod
    def _build_model_config(train_fdr, n_jobs, rng):
        return {"rescore": {"trainFDR": train_fdr, "numJobs": n_jobs, "seed": rng}}

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

        train_fdr = getattr(self, "train_fdr", 0.01)
        import optimhc.rescore.model  # noqa: F401 -- register configured models lazily
        from optimhc.rescore import mokapot as mokapot_adapter

        model_cls = rescore_model_factory.get_model(model_type)
        model = model_cls.from_config(
            self._build_model_config(train_fdr=train_fdr, n_jobs=n_jobs, rng=self.rng)
        )

        kwargs = {}
        if rescoring_features is not None:
            kwargs["rescoring_features"] = rescoring_features

        results, models = mokapot_adapter.rescore(
            psms,
            model=model,
            test_fdr=test_fdr,
            rng=self.rng,
            **kwargs,
        )
        return results, models

    def save_results(
        self,
        psms,
        results,
        models,
        output_dir=None,
        file_root="optimhc",
        feature_columns=None,
    ):
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
        feature_columns : sequence of str, optional
            Exact feature subset used for both rescoring and PIN serialization.
        """
        output_dir = output_dir if output_dir is not None else self.output_dir
        from optimhc.rescore import mokapot as mokapot_adapter

        results.to_txt(dest_dir=output_dir, file_root=file_root, decoys=True)
        mokapot_adapter.write_pin(
            psms,
            os.path.join(output_dir, f"{file_root}.pin"),
            feature_columns=feature_columns,
        )

        if self.save_models:
            model_dir = os.path.join(output_dir, "models")
            os.makedirs(model_dir, exist_ok=True)
            logger.info(f"Saving models to {model_dir}")
            for i, model in enumerate(models):
                model.save(os.path.join(model_dir, f"{file_root}.model{i}"))

        if self.to_flashlfq:
            from optimhc.output.flashlfq import write_flashlfq

            write_flashlfq(
                psms,
                results,
                os.path.join(output_dir, f"{file_root}.FlashLFQ.txt"),
                fdr=self.test_fdr,
            )

    def visualize_results(self, psms, results, models, output_dir=None):
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
        """
        if not self.visualization_enabled:
            logger.info("Visualization is disabled. Skipping...")
            return
        from optimhc.visualization import (
            plot_feature_importance,
            plot_qvalues,
            visualize_feature_correlation,
            visualize_target_decoy_features,
        )

        output_dir = output_dir if output_dir is not None else self.output_dir
        fig_dir = os.path.join(output_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)

        plot_qvalues(
            results,
            save_path=os.path.join(fig_dir, "qvalues.png"),
            threshold=0.05,
        )

        from optimhc.rescore.mokapot import mokapot_feature_columns

        plot_feature_importance(
            models,
            mokapot_feature_columns(psms),
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
        results = None
        models = None
        try:
            from optimhc.visualization import plot_feature_importance, plot_qvalues

            os.makedirs(exp_dir, exist_ok=True)

            model_type = exp_config.get("model", self.model_type)
            n_jobs = exp_config.get("numJobs", self.n_jobs)
            if "source" in exp_config:
                features = select_feature_groups(self.feature_groups, exp_config["source"])
            else:
                features = tuple(exp_config.get("features", psms.feature_columns))
            logger.info(f"Features used in experiment '{exp_name}': {features}")

            results, models = self.rescore(
                psms,
                model_type=model_type,
                n_jobs=n_jobs,
                test_fdr=self.test_fdr,
                rescoring_features=features,
            )

            self.save_results(
                psms,
                results,
                models,
                output_dir=exp_dir,
                file_root=exp_name,
                feature_columns=features,
            )

            fig_dir = os.path.join(exp_dir, "figures")

            plot_qvalues(
                results,
                save_path=os.path.join(fig_dir, "qvalues.png"),
                threshold=0.05,
            )

            from optimhc.rescore.mokapot import mokapot_feature_columns

            plot_feature_importance(
                models,
                feature_columns=mokapot_feature_columns(psms, features),
                save_path=os.path.join(fig_dir, "feature_importance.png"),
            )

            return True

        except Exception as e:
            logger.error(f"Experiment '{exp_name}' failed: {e}")
            return False

        finally:
            del results
            del models
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
        pin_path = os.path.join(self.output_dir, f"optimhc.{self.experiment}.pin")
        from optimhc.rescore import mokapot as mokapot_adapter
        from optimhc.visualization import visualize_feature_correlation

        mokapot_adapter.write_pin(psms, pin_path)
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
