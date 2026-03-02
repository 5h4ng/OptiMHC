"""
Pipeline control component for optiMHC GUI.
"""

import os
import subprocess
import sys
import tempfile
import time
from typing import Dict, Any, Optional, Tuple
import streamlit as st
import yaml

from optimhc.gui.utils import create_temp_config_file, run_pipeline_command


def pipeline_status_indicator(running: bool = False, success: Optional[bool] = None):
    """
    Display a status indicator for the pipeline.

    Args:
        running: Whether the pipeline is currently running
        success: Whether the pipeline completed successfully
    """
    if running:
        st.info("Pipeline is running...")
    elif success is not None:
        if success:
            st.success("Pipeline completed successfully")
        else:
            st.error("Pipeline failed")
    else:
        st.info("Pipeline not yet started")


def pipeline_control_panel(config: Dict[str, Any]):
    """
    Create a control panel for running the pipeline.

    Args:
        config: Configuration dictionary
    """
    st.subheader("Pipeline Control")

    # Initialize session state
    if "pipeline_running" not in st.session_state:
        st.session_state.pipeline_running = False

    if "pipeline_process" not in st.session_state:
        st.session_state.pipeline_process = None

    if "pipeline_start_time" not in st.session_state:
        st.session_state.pipeline_start_time = None

    if "pipeline_config_path" not in st.session_state:
        st.session_state.pipeline_config_path = None

    # Display status
    col1, col2 = st.columns([1, 3])

    with col1:
        if st.session_state.pipeline_running:
            pipeline_status_indicator(running=True)
        else:
            if st.session_state.pipeline_process is not None:
                return_code = st.session_state.pipeline_process.poll()
                pipeline_status_indicator(success=(return_code == 0))
            else:
                pipeline_status_indicator()

    with col2:
        if st.session_state.pipeline_start_time:
            elapsed_time = time.time() - st.session_state.pipeline_start_time
            st.text(f"Running for: {int(elapsed_time // 60)}m {int(elapsed_time % 60)}s")

    # Control buttons
    start_disabled = st.session_state.pipeline_running
    stop_disabled = not st.session_state.pipeline_running

    col1, col2 = st.columns(2)

    with col1:
        if st.button("Start Pipeline", disabled=start_disabled, key="start_pipeline"):
            # Check if the configuration is valid
            if not config.get("experimentName"):
                st.error("Experiment name is required")
                return

            if not config.get("inputFile"):
                st.error("At least one input file is required")
                return

            if not config.get("allele"):
                st.error("At least one allele is required")
                return

            if not config.get("featureGenerator"):
                st.error("At least one feature generator is required")
                return

            # Create a temporary configuration file
            config_path = create_temp_config_file(config)
            st.session_state.pipeline_config_path = config_path

            # Run the pipeline as a subprocess
            st.session_state.pipeline_process = run_pipeline_command(config_path)
            st.session_state.pipeline_running = True
            st.session_state.pipeline_start_time = time.time()

            # Initialize logs
            if "logs" not in st.session_state:
                st.session_state.logs = []

            # Rerun to update UI
            st.rerun()

    with col2:
        if st.button("Stop Pipeline", disabled=stop_disabled, key="stop_pipeline"):
            if st.session_state.pipeline_process:
                # Terminate the process
                st.session_state.pipeline_process.terminate()
                st.session_state.pipeline_running = False

                # Wait for process to terminate
                try:
                    st.session_state.pipeline_process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    # Force kill if it doesn't terminate gracefully
                    st.session_state.pipeline_process.kill()

                st.warning("Pipeline was stopped by user")

                # No longer need to cleanup config file since it's part of the output

                # Rerun to update UI
                st.rerun()

    # Save configuration button
    if st.button("Save Configuration to File"):
        # Create a download button for the configuration
        config_yaml = yaml.dump(config, default_flow_style=False)

        # Use streamlit's download button
        filename = f"{config.get('experimentName', 'optimhc_config')}.yaml"
        st.download_button(
            label="Download Configuration File",
            data=config_yaml,
            file_name=filename,
            mime="text/yaml"
        )


def check_pipeline_status():
    """
    Check the status of a running pipeline.

    Returns:
        Tuple of (running, return_code)
    """
    running = st.session_state.get("pipeline_running", False)
    process = st.session_state.get("pipeline_process", None)

    if process is None:
        return False, None

    # Check if process is still running
    return_code = process.poll()

    if return_code is not None and running:
        # Process has completed
        st.session_state.pipeline_running = False

        # No longer need to cleanup config file since it's part of the output

        return False, return_code

    return running, return_code
