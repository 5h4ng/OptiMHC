"""
Utility functions for the optiMHC GUI.
"""

import os
import subprocess
import sys
import tempfile
from pathlib import Path
import yaml
import json
import streamlit as st
import pandas as pd
from typing import Dict, List, Any, Optional, Union


def load_config_from_yaml(file_path: str) -> Dict[str, Any]:
    """
    Load configuration from a YAML file.

    Args:
        file_path: Path to the YAML configuration file

    Returns:
        Dictionary containing the configuration
    """
    try:
        with open(file_path, 'r') as f:
            config = yaml.safe_load(f)
        return config
    except Exception as e:
        st.error(f"Error loading configuration file: {str(e)}")
        return {}


def save_config_to_yaml(config: Dict[str, Any], file_path: str) -> bool:
    """
    Save configuration to a YAML file.

    Args:
        config: Configuration dictionary
        file_path: Path to save the YAML file

    Returns:
        True if successful, False otherwise
    """
    try:
        with open(file_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False)
        return True
    except Exception as e:
        st.error(f"Error saving configuration file: {str(e)}")
        return False


def create_temp_config_file(config: Dict[str, Any]) -> str:
    """
    Create a temporary configuration file.

    Args:
        config: Configuration dictionary

    Returns:
        Path to the temporary file
    """
    fd, temp_path = tempfile.mkstemp(suffix='.yaml')
    with os.fdopen(fd, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    return temp_path


def run_pipeline_command(config_path: str) -> subprocess.Popen:
    """
    Run the optiMHC pipeline as a subprocess.

    Args:
        config_path: Path to the configuration file

    Returns:
        Subprocess object
    """
    # Load configuration to access output directory and experiment name
    config = {}
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        # Store config in session state for log finding
        st.session_state.config = config
    except Exception as e:
        print(f"Error loading configuration: {str(e)}")

    # Identify where logs will be stored based on pipeline conventions
    output_dir = config.get("outputDir", "")
    experiment_name = config.get("experimentName", "")

    # The actual log file that the pipeline will create
    if output_dir and experiment_name:
        experiment_dir = os.path.join(output_dir, experiment_name)
        os.makedirs(experiment_dir, exist_ok=True)
        pipeline_log_path = os.path.join(experiment_dir, "log")

        # Store the expected log file path in session state
        st.session_state.pipeline_log_path = pipeline_log_path
    else:
        st.session_state.pipeline_log_path = None

    # Reset log position counter
    st.session_state.log_position = 0

    # Set environment variables to ensure output is not buffered
    my_env = os.environ.copy()
    my_env["PYTHONUNBUFFERED"] = "1"  # Prevent Python from buffering output
    my_env["PYTHONIOENCODING"] = "utf-8"  # Ensure UTF-8 encoding

    # Run the pipeline command
    cmd = [sys.executable, "-u", "-m", "optimhc", "pipeline", "--config", config_path]
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,  # Line buffered
        universal_newlines=True,
        env=my_env
    )

    # Start a thread to read and print output (for debugging)
    def monitor_output():
        try:
            for line in iter(process.stdout.readline, ''):
                print(line.strip())  # Print to console for debugging
        except Exception as e:
            print(f"Error in monitor thread: {str(e)}")

    # Start thread in daemon mode so it won't block program exit
    import threading
    monitor_thread = threading.Thread(target=monitor_output, daemon=True)
    monitor_thread.start()

    return process


def scan_output_directory(output_dir: str) -> Dict[str, List[str]]:
    """
    Scan the output directory for results.

    Args:
        output_dir: Path to the output directory

    Returns:
        Dictionary with lists of files grouped by type
    """
    if not os.path.exists(output_dir):
        return {
            'figures': [],
            'tables': [],
            'logs': [],
            'other': []
        }

    result_files = {
        'figures': [],
        'tables': [],
        'logs': [],
        'other': []
    }

    for root, _, files in os.walk(output_dir):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith(('.png', '.jpg', '.jpeg', '.svg', '.pdf')):
                result_files['figures'].append(file_path)
            elif file.endswith(('.csv', '.tsv', '.xlsx')):
                result_files['tables'].append(file_path)
            elif file.endswith('.log'):
                result_files['logs'].append(file_path)
            else:
                result_files['other'].append(file_path)

    return result_files


def get_config_summary(config: Dict[str, Any]) -> Dict[str, Any]:
    """
    Create a summary of the configuration for display.

    Args:
        config: Configuration dictionary

    Returns:
        Dictionary with summarized configuration
    """
    summary = {
        'Experiment Name': config.get('experimentName', 'N/A'),
        'Input Type': config.get('inputType', 'N/A'),
        'Output Directory': config.get('outputDir', 'N/A'),
        'Alleles': ', '.join(config.get('allele', [])) if config.get('allele') else 'N/A',
        'Feature Generators': [gen.get('name') for gen in config.get('featureGenerator', [])]
                              if config.get('featureGenerator') else 'N/A',
        'Rescoring Model': config.get('rescore', {}).get('model', 'N/A') if config.get('rescore') else 'N/A',
        'Modification Map': f"{len(config.get('modificationMap', {}))} modifications"
                          if config.get('modificationMap') else 'Default'
    }
    return summary


def parse_feature_generator_json(json_str: str) -> Optional[Dict[str, Any]]:
    """
    Parse JSON string for feature generator configuration.

    Args:
        json_str: JSON string to parse

    Returns:
        Dictionary containing the feature generator configuration or None if invalid
    """
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        return None


def format_log_message(message: str, level: str = "INFO") -> str:
    """
    Format a log message for display.

    Args:
        message: Log message to format
        level: Log level

    Returns:
        Formatted log message
    """
    level_colors = {
        "DEBUG": "gray",
        "INFO": "black",
        "WARNING": "orange",
        "ERROR": "red",
        "CRITICAL": "crimson"
    }
    color = level_colors.get(level, "black")
    return f'<span style="color: {color}">[{level}] {message}</span>'
