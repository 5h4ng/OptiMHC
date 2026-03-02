"""
Results page for the optiMHC GUI.
"""

import streamlit as st
import os
import glob
from typing import Dict, Any, List

from optimhc.gui.style import main_header, sub_header, info_box, warning_box
from optimhc.gui.components.results_viewer import results_viewer


def find_figure_files(output_dir: str, experiment_name: str = None) -> List[str]:
    """
    Find all figure files in the output directory or its figures subdirectory.

    Args:
        output_dir: Path to the output directory
        experiment_name: Optional experiment name subdirectory

    Returns:
        List of paths to figure files
    """
    figure_files = []

    # If experiment name is provided, look in that subdirectory
    if experiment_name:
        experiment_dir = os.path.join(output_dir, experiment_name)
        if os.path.exists(experiment_dir):
            # Look for figures subfolder inside experiment directory
            figures_dir = os.path.join(experiment_dir, "figures")
            if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
                figure_files = glob.glob(os.path.join(figures_dir, "*.png"))
                figure_files.extend(glob.glob(os.path.join(figures_dir, "*.jpg")))
                figure_files.extend(glob.glob(os.path.join(figures_dir, "*.svg")))
                figure_files.extend(glob.glob(os.path.join(figures_dir, "*.pdf")))

            # If no figures found in subfolder, look in experiment directory
            if not figure_files:
                figure_files = glob.glob(os.path.join(experiment_dir, "*.png"))
                figure_files.extend(glob.glob(os.path.join(experiment_dir, "*.jpg")))
                figure_files.extend(glob.glob(os.path.join(experiment_dir, "*.svg")))
                figure_files.extend(glob.glob(os.path.join(experiment_dir, "*.pdf")))

    # If no figures found in experiment directory or no experiment name,
    # check main output directory
    if not figure_files:
        # Check for figures subfolder in main directory
        figures_dir = os.path.join(output_dir, "figures")
        if os.path.exists(figures_dir) and os.path.isdir(figures_dir):
            figure_files = glob.glob(os.path.join(figures_dir, "*.png"))
            figure_files.extend(glob.glob(os.path.join(figures_dir, "*.jpg")))
            figure_files.extend(glob.glob(os.path.join(figures_dir, "*.svg")))
            figure_files.extend(glob.glob(os.path.join(figures_dir, "*.pdf")))

        # If still no figures, check main directory
        if not figure_files:
            figure_files = glob.glob(os.path.join(output_dir, "*.png"))
            figure_files.extend(glob.glob(os.path.join(output_dir, "*.jpg")))
            figure_files.extend(glob.glob(os.path.join(output_dir, "*.svg")))
            figure_files.extend(glob.glob(os.path.join(output_dir, "*.pdf")))

    return figure_files


def render():
    """
    Render the results page focusing on figure visualization.
    """
    main_header("Results Visualization")

    # Check if we have a configuration
    if "config" not in st.session_state or not st.session_state.config:
        warning_box("No configuration available. Please configure and run the pipeline first.")

        if st.button("Go to Configuration"):
            st.session_state.page = "configure"
            st.rerun()
        return

    # Get output directory and experiment name from configuration
    output_dir = st.session_state.config.get("outputDir", "./results")
    experiment_name = st.session_state.config.get("experimentName", "")

    # Allow user to change the output directory
    output_dir = st.text_input("Output Directory", value=output_dir)

    # Check if directory exists
    if not os.path.exists(output_dir):
        warning_box(f"Output directory '{output_dir}' does not exist. Please run the pipeline first or check the directory path.")

        # Show navigation buttons
        col1, col2 = st.columns(2)
        with col1:
            if st.button("← Back to Run Pipeline"):
                st.session_state.page = "run"
                st.rerun()
        with col2:
            if st.button("Go to Home"):
                st.session_state.page = "home"
                st.rerun()
        return

    # Show experiment details if available
    if experiment_name:
        sub_header(f"Experiment: {experiment_name}")

    # Display results using our results_viewer component
    # The component will handle finding the appropriate figure files
    results_viewer(output_dir)

    # Navigation buttons
    st.markdown("---")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("← Back to Run"):
            st.session_state.page = "run"
            st.rerun()

    with col2:
        if st.button("Back to Home"):
            st.session_state.page = "home"
            st.rerun()
