import logging
import os
from copy import deepcopy

import yaml

logger = logging.getLogger(__name__)

# Default configuration with all possible parameters
DEFAULT_CONFIG = {
    "outputDir": "./results",
    "inputType": "pepxml",
    "inputFile": [],
    "decoyPrefix": "DECOY_",
    "visualization": True,
    "saveModels": True,
    "toFlashLFQ": True,
    "allele": [],
    "numProcess": 4,
    "removePreNxtAA": False,
    "showProgress": True,
    "logLevel": "INFO",
    "rescore": {"testFDR": 0.01, "trainFDR": 0.01, "model": "Percolator", "numJobs": 1},
}


def _deep_merge(default, override):
    """
    Deep merge two dictionaries. The override dictionary values take precedence.

    Parameters
    ----------
    default : dict
        Default dictionary.
    override : dict
        Dictionary with values to override defaults.

    Returns
    -------
    dict
        Merged dictionary.
    """
    result = deepcopy(default)

    if not isinstance(override, dict):
        return override

    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value

    return result


def load_config(config_path):
    """
    Load and parse a configuration file using YAML.
    Merges loaded config with default configuration.

    Parameters:
        config_path (str): Path to the YAML configuration file.

    Returns:
        dict: A dictionary containing all configurations.
    """
    logger.info(f"Loading configuration from {config_path}")
    with open(config_path, "r") as f:
        user_config = yaml.safe_load(f)

    config = _deep_merge(DEFAULT_CONFIG, user_config)

    return config


class Config:
    """
    Configuration manager for optiMHC pipeline.

    This class handles loading, validating, and providing access to configuration parameters
    from YAML files or dictionaries. It implements a fail-fast validation strategy to ensure
    configuration correctness before pipeline execution.

    Parameters
    ----------
    config_source : str or dict or None, optional
        Path to YAML file, dictionary with configuration, or None for default config.
        If None, uses DEFAULT_CONFIG.

    Attributes
    ----------
    _config : dict
        The internal configuration dictionary.

    Raises
    ------
    ValueError
        If configuration is invalid or required parameters are missing.
    FileNotFoundError
        If specified YAML file does not exist.
    yaml.YAMLError
        If YAML file is malformed.

    Examples
    --------
    >>> # Load from YAML file
    >>> config = Config("path/to/config.yaml")

    >>> # Load from dictionary
    >>> config_dict = {
    ...     "inputType": "pepxml",
    ...     "inputFile": ["data.pep.xml"],
    ...     "outputDir": "./results",
    ...     "rescore": {"testFDR": 0.01, "model": "Percolator"}
    ... }
    >>> config = Config(config_dict)

    >>> # Use default configuration
    >>> config = Config()

    >>> # Access configuration values
    >>> input_type = config["inputType"]
    >>> output_dir = config.get("outputDir", "./default")

    >>> # Save configuration to file
    >>> config.save("output_config.yaml")
    """

    def __init__(self, config_source=None):
        """
        Initialize Config from a YAML file path or a dictionary.

        Parameters
        ----------
        config_source : str or dict or None, optional
            Path to YAML file, dictionary with configuration, or None for default config.
            If None, uses DEFAULT_CONFIG.

        Raises
        ------
        ValueError
            If config_source is neither a string, dict, nor None.
        FileNotFoundError
            If specified YAML file does not exist.
        yaml.YAMLError
            If YAML file is malformed.
        """
        if config_source is None:
            self._config = deepcopy(DEFAULT_CONFIG)
        elif isinstance(config_source, str):
            with open(config_source, "r") as f:
                user_config = yaml.safe_load(f)
            self._config = _deep_merge(DEFAULT_CONFIG, user_config)
        elif isinstance(config_source, dict):
            self._config = _deep_merge(DEFAULT_CONFIG, config_source)
        else:
            raise ValueError("Config source must be a file path, dict, or None.")

    # TODO: allele name validation
    def validate(self):
        """
        Validate the configuration using a fail-fast strategy.

        This method performs comprehensive validation of the configuration,
        including required fields, data types, file existence, and feature
        generator configurations.

        Raises
        ------
        ValueError
            If any validation check fails. The error message will indicate
            the specific validation failure.

        Notes
        -----
        The validation includes checks for:
        - Required fields (inputType, inputFile, outputDir, rescore)
        - Input file existence and type
        - Output directory creation
        - Rescore configuration validity (TODO)
        - Feature generator configuration primitive validity (TODO)
        - Optional parameter validity (TODO): we should validate 'allele' first !!!
        """
        if not isinstance(self._config, dict):
            logger.error("Configuration must be a dictionary")
            raise ValueError("Configuration must be a dictionary")

        required_fields = ["inputType", "inputFile", "outputDir", "rescore"]
        for field in required_fields:
            if field not in self._config:
                logger.error(f"Missing required configuration: '{field}'")
                raise ValueError(f"Missing required configuration: '{field}'")

            if field == "inputFile" and self._config[field] == []:
                logger.error("inputFile list cannot be empty")
                raise ValueError("inputFile list cannot be empty")
            elif self._config[field] in (None, "", []):
                logger.error(f"Required configuration '{field}' cannot be empty")
                raise ValueError(f"Required configuration '{field}' cannot be empty")

        if self._config["inputType"] not in ("pepxml", "pin"):
            logger.error("inputType must be 'pepxml' or 'pin'")
            raise ValueError("inputType must be 'pepxml' or 'pin'")

        if (
            self._config["inputType"] == "pin"
            and self._config.get("retentionTimeColumn", None) is None
        ):
            logger.error("retentionTimeColumn must be specified when inputType is 'pin'")
            raise ValueError("retentionTimeColumn must be specified when inputType is 'pin'")

        input_files = self._config["inputFile"]
        if not isinstance(input_files, (list, tuple)):
            logger.debug(f"inputFile is not a list or tuple: {input_files}. Converting to list.")
            self._config["inputFile"] = list(input_files)
        if not input_files:
            logger.error("inputFile list cannot be empty")
            raise ValueError("inputFile list cannot be empty")

        for file_path in input_files:
            if not os.path.exists(file_path):
                logger.error(f"Input file does not exist: {file_path}")
                raise ValueError(f"Input file does not exist: {file_path}")

        output_dir = self._config["outputDir"]
        if not isinstance(output_dir, str):
            logger.error("outputDir must be a string")
            raise ValueError("outputDir must be a string")
        if not output_dir:
            logger.error("outputDir is required")
            raise ValueError("outputDir is required")
        os.makedirs(output_dir, exist_ok=True)

        # TODO: Validate feature generator configuration
        valid_generators = {
            "Basic",
            "OverlappingPeptide",
            "PWM",
            "MHCflurry",
            "NetMHCpan",
            "NetMHCIIpan",
            "DeepLC",
            "SpectralSimilarity",
        }

        if self._config["featureGenerator"] is None:
            logger.warning("No feature generators specified.")
            self._config["featureGenerator"] = []
        for fg in self._config["featureGenerator"]:
            if fg["name"] not in valid_generators:
                logger.error(f"Invalid feature generator: {fg['name']}")
                raise ValueError(f"Invalid feature generator: {fg['name']}")

        valid_sources = valid_generators | {"Original", "ContigFeatures"}
        if self._config.get("experiments", None) is not None:  # experiment mode
            if not isinstance(self._config["experiments"], list):
                logger.error("experiments must be a list")
                raise ValueError("experiments must be a list")
            for exp in self._config["experiments"]:
                for source in exp.get("source", []):
                    if source not in valid_sources:
                        logger.error(f"Invalid source in experiments: {source}")
                        raise ValueError(f"Invalid source in experiments: {source}")

    def to_dict(self):
        return deepcopy(self._config)

    def save(self, path):
        with open(path, "w") as f:
            yaml.safe_dump(self._config, f)

    def __getitem__(self, key):
        return self._config[key]

    def __setitem__(self, key, value):
        self._config[key] = value

    def get(self, key, default=None):
        return self._config.get(key, default)

    def __contains__(self, key):
        return key in self._config

    def __repr__(self):
        return f"Config({self._config})"
