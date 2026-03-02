"""
Run page for the optiMHC GUI.
"""

import streamlit as st
import os
import time
from typing import Dict, Any

from optimhc.gui.style import main_header, sub_header, info_box, warning_box
from optimhc.gui.components.pipeline_control import pipeline_control_panel, check_pipeline_status
from optimhc.gui.components.log_viewer import log_viewer
from optimhc.gui.components.config_form import render_config_summary


def render():
    """
    Render the run page.
    """
    main_header("Run Pipeline")

    # Check if we have a configuration
    if "config" not in st.session_state or not st.session_state.config:
        warning_box("No configuration available. Please configure the pipeline first.")

        if st.button("Go to Configuration"):
            st.session_state.page = "configure"
            st.rerun()
        return

    # Display configuration summary
    with st.expander("Configuration Summary", expanded=False):
        render_config_summary(st.session_state.config)

    # Pipeline control panel
    pipeline_control_panel(st.session_state.config)

    # Separator
    st.markdown("---")

    # Check pipeline status
    running, return_code = check_pipeline_status()

    # Show log viewer if process exists
    if "pipeline_process" in st.session_state and st.session_state.pipeline_process:
        log_viewer(st.session_state.pipeline_process)

    # Navigation buttons
    col1, col2 = st.columns(2)

    with col1:
        if st.button("← Back to Configure"):
            st.session_state.page = "configure"
            st.rerun()

    with col2:
        # Only enable results button if pipeline has completed
        results_disabled = running or ("pipeline_process" not in st.session_state)

        if st.button("View Results →", disabled=results_disabled):
            st.session_state.page = "results"
            st.rerun()
