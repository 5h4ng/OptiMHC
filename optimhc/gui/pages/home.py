"""
Home page for the optiMHC GUI.
"""

import streamlit as st
from optimhc.gui.style import main_header, sub_header, info_box


def render():
    """
    Render the home page.
    """
    main_header("Welcome to optiMHC")

    st.markdown("""
    optiMHC is a high-performance rescoring pipeline for immunopeptidomics data, designed to enhance
    peptide identification by integrating multiple feature generators and machine learning-based rescoring.
    """)

    # Feature overview
    sub_header("Features")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        - **Multiple Feature Generators**
          - Basic sequence features
          - Position Weight Matrix (PWM)
          - MHC binding prediction (MHCflurry, NetMHCpan)
          - Chromatographic retention time (DeepLC)
          - Spectrum similarity
          - Overlapping peptides

        - **Machine Learning Rescoring**
          - Percolator
          - XGBoost
          - Random Forest
        """)

    with col2:
        st.markdown("""
        - **Flexible Configuration**
          - YAML-based configuration
          - Command-line interface
          - Python API

        - **Comprehensive Visualization**
          - ROC curves
          - Score distributions
          - FDR vs. identifications
          - Feature importance
        """)

    # Quick start guide
    sub_header("Quick Start")

    st.markdown("""
    To get started with optiMHC:

    1. Go to the **Configure** page to set up your pipeline parameters
    2. Upload your input files or specify their paths
    3. Run the pipeline from the **Run** page
    4. View your results in the **Results** page
    """)

    # Warning about GUI being in development
    info_box("""
    <strong>Note:</strong> This GUI is currently in development. For advanced use cases or troubleshooting,
    please refer to the command-line interface or the documentation.
    """)

    # Footer
    st.markdown("---")
    st.markdown("optiMHC is developed and maintained by Zixiang Shang. For more information, please visit the [GitHub repository](https://github.com/5h4ng/optiMHC).")
