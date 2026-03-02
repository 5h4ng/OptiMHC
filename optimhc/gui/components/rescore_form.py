"""
Rescoring form component for optiMHC GUI.
"""

import streamlit as st
from typing import Dict, Any

# Import optiMHC config defaults
from optimhc.core.config import DEFAULT_CONFIG

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
