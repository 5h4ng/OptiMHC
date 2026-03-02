"""
Main application for the optiMHC GUI.
"""

import streamlit as st
import os
import sys
from pathlib import Path

# Add optiMHC root to path if needed
optimhc_root = str(Path(__file__).parent.parent.parent)
if optimhc_root not in sys.path:
    sys.path.append(optimhc_root)

# Import style utilities
from optimhc.gui.style import set_page_config, apply_custom_css, footer

# Import page modules
from optimhc.gui.pages import home, configure, run, results


def main():
    """
    Main application entry point.
    """
    # Set up page config
    set_page_config()

    # Apply custom CSS
    apply_custom_css()

    # Initialize session state for navigation
    if "page" not in st.session_state:
        st.session_state.page = "home"

    # Sidebar navigation
    st.sidebar.title("Navigation")

    # Navigation buttons
    if st.sidebar.button("Home", use_container_width=True):
        st.session_state.page = "home"
        st.rerun()

    if st.sidebar.button("Configure", use_container_width=True):
        st.session_state.page = "configure"
        st.rerun()

    if st.sidebar.button("Run Pipeline", use_container_width=True):
        st.session_state.page = "run"
        st.rerun()

    if st.sidebar.button("Results", use_container_width=True):
        st.session_state.page = "results"
        st.rerun()

    # Version info in sidebar
    st.sidebar.markdown("---")
    from optimhc import __version__
    st.sidebar.info(f"optiMHC v{__version__}")

    # Render the selected page
    if st.session_state.page == "home":
        home.render()
    elif st.session_state.page == "configure":
        configure.render()
    elif st.session_state.page == "run":
        run.render()
    elif st.session_state.page == "results":
        results.render()

    # Footer
    footer()


if __name__ == "__main__":
    main()
