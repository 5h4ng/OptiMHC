"""
Configuration form component for optiMHC GUI.
"""

import os
import streamlit as st
import yaml
from typing import Dict, Any, List, Optional
import json

# Import optiMHC config defaults
from optimhc.core.config import DEFAULT_CONFIG


def feature_generator_form(existing_generators: List[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """
    Create a form for configuring feature generators.

    Args:
        existing_generators: List of existing feature generator configurations

    Returns:
        List of feature generator configurations
    """
    if existing_generators is None:
        existing_generators = []

    feature_generators = []

    # Known feature generators and their parameters
    generator_options = {
        "Basic": {},
        "PWM": {"class": ["I", "II"]},
        "MHCflurry": {},
        "NetMHCpan": {},
        "NetMHCIIpan": {},
        "DeepLC": {
            "calibrationCriteria": ["expect", "xcorr", "hyperscore"],
            "lowerIsBetter": [True, False],
            "calibrationSize": [0.1, 0.2, 0.3]
        },
        "SpectralSimilarity": {
            "model": ["AlphaPeptDeep_ms2_generic", "AlphaPeptDeep_ms2_HCD", "AlphaPeptDeep_ms2_CID"],
            "instrument": ["LUMOS", "QE", "VELOS", "FUSION"],
            "numTopPeaks": [10, 20, 36, 50]
        },
        "OverlappingPeptide": {
            "minOverlapLength": [7, 8, 9],
            "minLength": [7, 8, 9],
            "maxLength": [15, 20, 25],
            "overlappingScore": ["expect", "xcorr", "hyperscore"]
        }
    }

    st.subheader("Feature Generators")

    # Use session state to keep track of the number of generators
    if "num_generators" not in st.session_state:
        st.session_state.num_generators = max(1, len(existing_generators))

    # Add/remove generator controls outside of the form
    col1, col2 = st.columns([1, 5])
    with col1:
        if st.button("➕ Add Generator", key="add_generator"):
            st.session_state.num_generators += 1
            st.rerun()
    with col2:
        if st.session_state.num_generators > 1 and st.button("➖ Remove Last Generator", key="remove_generator"):
            st.session_state.num_generators -= 1
            st.rerun()

    # Generate forms for each feature generator
    for i in range(st.session_state.num_generators):
        with st.expander(f"Feature Generator {i+1}", expanded=True):
            existing_gen = {} if i >= len(existing_generators) else existing_generators[i]

            # Feature generator name
            generator_name = st.selectbox(
                "Generator Type",
                options=list(generator_options.keys()),
                key=f"gen_type_{i}",
                index=list(generator_options.keys()).index(existing_gen.get("name", "Basic")) if existing_gen.get("name") in generator_options else 0
            )

            # Feature generator parameters
            params = {}
            if generator_options[generator_name]:
                st.markdown("**Parameters:**")
                for param_name, param_options in generator_options[generator_name].items():
                    existing_params = existing_gen.get("params", {})
                    existing_value = existing_params.get(param_name)

                    # Handle different parameter types
                    if isinstance(param_options, list):
                        if all(isinstance(x, bool) for x in param_options):
                            param_value = st.checkbox(
                                param_name,
                                value=existing_value if existing_value is not None else param_options[0],
                                key=f"gen_{i}_{param_name}"
                            )
                        elif all(isinstance(x, (int, float)) for x in param_options):
                            param_value = st.number_input(
                                param_name,
                                value=existing_value if existing_value is not None else param_options[0],
                                key=f"gen_{i}_{param_name}"
                            )
                        else:
                            param_value = st.selectbox(
                                param_name,
                                options=param_options,
                                index=param_options.index(existing_value) if existing_value in param_options else 0,
                                key=f"gen_{i}_{param_name}"
                            )
                    else:
                        param_value = st.text_input(
                            param_name,
                            value=str(existing_value) if existing_value is not None else "",
                            key=f"gen_{i}_{param_name}"
                        )

                    params[param_name] = param_value

            # Add to list of generators
            generator_config = {"name": generator_name}
            if params:
                generator_config["params"] = params

            feature_generators.append(generator_config)

    return feature_generators


def rescore_form(existing_rescore: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Create a form for rescoring settings.

    Args:
        existing_rescore: Existing rescore configuration

    Returns:
        Rescore configuration dictionary
    """
    if existing_rescore is None:
        existing_rescore = DEFAULT_CONFIG["rescore"]

    st.subheader("Rescoring Settings")

    rescore_model = st.selectbox(
        "Rescoring Model",
        options=["Percolator", "XGBoost", "RandomForest"],
        index=["Percolator", "XGBoost", "RandomForest"].index(existing_rescore.get("model", "Percolator")),
        help="Model to use for rescoring"
    )

    test_fdr = st.number_input(
        "Test FDR",
        min_value=0.001,
        max_value=0.1,
        value=float(existing_rescore.get("testFDR", 0.01)),
        step=0.001,
        format="%.3f",
        help="FDR threshold for testing"
    )

    num_jobs = st.number_input(
        "Number of Jobs",
        min_value=1,
        max_value=32,
        value=int(existing_rescore.get("numJobs", 1)),
        help="Number of parallel jobs for model training"
    )

    return {
        "model": rescore_model,
        "testFDR": test_fdr,
        "numJobs": num_jobs
    }


def config_form(existing_config: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Create a form for configuring the pipeline.

    Args:
        existing_config: Existing configuration dictionary

    Returns:
        Configuration dictionary
    """
    if existing_config is None:
        existing_config = DEFAULT_CONFIG

    st.subheader("Basic Settings")

    experiment_name = st.text_input(
        "Experiment Name",
        value=existing_config.get("experimentName", ""),
        help="Name of the experiment"
    )

    input_type = st.selectbox(
        "Input Type",
        options=["pepxml", "pin"],
        index=["pepxml", "pin"].index(existing_config.get("inputType", "pepxml")),
        help="Type of input file"
    )

    # For GUI, we'll handle input files differently than the direct file paths
    input_files = existing_config.get("inputFile", [])
    if isinstance(input_files, str):
        input_files = [input_files]

    input_files_str = st.text_area(
        "Input Files",
        value="\n".join(input_files) if input_files else "",
        height=100,
        help="One file path per line. Use file uploader to add files."
    )

    input_files = [f for f in input_files_str.strip().split("\n") if f]

    decoy_prefix = st.text_input(
        "Decoy Prefix",
        value=existing_config.get("decoyPrefix", "DECOY_"),
        help="Prefix used to identify decoy sequences"
    )

    output_dir = st.text_input(
        "Output Directory",
        value=existing_config.get("outputDir", "./results"),
        help="Directory where results will be saved"
    )

    # Allele settings
    st.subheader("Allele Settings")

    alleles = existing_config.get("allele", [])
    if isinstance(alleles, str):
        alleles = [alleles]

    alleles_str = st.text_area(
        "Alleles",
        value="\n".join(alleles) if alleles else "",
        height=100,
        help="One allele per line, e.g., HLA-A*02:01"
    )

    alleles = [a for a in alleles_str.strip().split("\n") if a]

    # Performance settings
    st.subheader("Performance Settings")

    col1, col2 = st.columns(2)

    with col1:
        num_processes = st.number_input(
            "Number of Processes",
            min_value=1,
            max_value=64,
            value=int(existing_config.get("numProcesses", 4)),
            help="Number of parallel processes"
        )

    with col2:
        show_progress = st.checkbox(
            "Show Progress",
            value=existing_config.get("showProgress", True),
            help="Show progress bars during processing"
        )

    col1, col2 = st.columns(2)

    with col1:
        visualization = st.checkbox(
            "Enable Visualization",
            value=existing_config.get("visualization", True),
            help="Generate visualizations of results"
        )

    with col2:
        remove_pre_nxt_aa = st.checkbox(
            "Remove Pre/Next Amino Acids",
            value=existing_config.get("removePreNxtAA", False),
            help="Remove pre/post neighboring amino acids in sequence processing"
        )

    log_level = st.selectbox(
        "Log Level",
        options=["DEBUG", "INFO", "WARNING", "ERROR"],
        index=["DEBUG", "INFO", "WARNING", "ERROR"].index(existing_config.get("logLevel", "INFO")),
        help="Logging verbosity level"
    )

    # Advanced sections

    # Feature generators
    feature_generators = feature_generator_form(existing_config.get("featureGenerator", []))

    # Rescoring
    rescore = rescore_form(existing_config.get("rescore", {}))

    # Combine all settings
    config = {
        "experimentName": experiment_name,
        "inputType": input_type,
        "inputFile": input_files,
        "decoyPrefix": decoy_prefix,
        "outputDir": output_dir,
        "allele": alleles,
        "numProcesses": num_processes,
        "showProgress": show_progress,
        "visualization": visualization,
        "removePreNxtAA": remove_pre_nxt_aa,
        "logLevel": log_level,
        "featureGenerator": feature_generators,
        "rescore": rescore
    }

    return config


def render_config_summary(config: Dict[str, Any]):
    """
    Render a summary of the configuration as a YAML code block.

    Args:
        config: Configuration dictionary
    """
    st.subheader("Configuration Summary")

    # Create a simplified copy of the configuration to display
    display_config = config.copy()

    # Convert to YAML string
    config_yaml = yaml.dump(display_config, default_flow_style=False, sort_keys=False)

    # Display as a code block with syntax highlighting
    st.code(config_yaml, language="yaml")


def validate_config(config: Dict[str, Any]) -> List[str]:
    """
    Validate configuration for obvious errors.

    Args:
        config: Configuration dictionary

    Returns:
        List of error messages, empty if configuration is valid
    """
    errors = []

    # Check required fields
    required_fields = ["experimentName", "inputType", "inputFile", "outputDir", "allele"]
    for field in required_fields:
        if field not in config or not config[field]:
            errors.append(f"Missing required field: {field}")

    # Check inputType
    if config.get("inputType") not in ["pepxml", "pin"]:
        errors.append("Input type must be 'pepxml' or 'pin'")

    # Check feature generators
    generators = config.get("featureGenerator", [])
    if not generators:
        errors.append("At least one feature generator is required")

    # Check if SpectralSimilarity has required parameters
    for gen in generators:
        if gen.get("name") == "SpectralSimilarity":
            params = gen.get("params", {})
            # Check instrument
            instrument = params.get("instrument")
            if instrument and instrument not in ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]:
                errors.append(f"Invalid instrument '{instrument}' for SpectralSimilarity. Must be one of: QE, LUMOS, TIMSTOF, SCIEXTOF")

            # Check mzML directory
            if "mzmlDir" not in params:
                errors.append("SpectralSimilarity requires 'mzmlDir' parameter")

            # Check spectrum ID pattern
            if "spectrumIdPattern" not in params:
                errors.append("SpectralSimilarity requires 'spectrumIdPattern' parameter to extract mzML filenames from spectrum IDs")

    # Check rescore settings
    rescore = config.get("rescore", {})
    if not rescore or "model" not in rescore:
        errors.append("Rescore model is required")

    if "testFDR" in rescore and (rescore["testFDR"] <= 0 or rescore["testFDR"] > 1):
        errors.append("Test FDR must be between 0 and 1")

    # Check modification map format
    mod_map = config.get("modificationMap", {})
    for mass, unimod in mod_map.items():
        if not unimod.startswith("UNIMOD:"):
            errors.append(f"Invalid UNIMOD format for mass {mass}: {unimod}. Must start with 'UNIMOD:'")

    return errors
