"""
Styling utilities for the optiMHC Streamlit interface.
"""

import streamlit as st


def apply_custom_css():
    """
    Apply custom CSS styling to the Streamlit app.
    """
    st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        color: #4257B2;
        margin-bottom: 1rem;
    }

    .sub-header {
        font-size: 1.5rem;
        color: #5273E0;
        margin-bottom: 0.5rem;
    }

    .info-box {
        background-color: #f0f2f6;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
    }

    .success-box {
        background-color: #d1f0d5;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
    }

    .warning-box {
        background-color: #fff4e0;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
    }

    .error-box {
        background-color: #ffe0e0;
        border-radius: 0.5rem;
        padding: 1rem;
        margin-bottom: 1rem;
    }

    .log-container {
        background-color: #f0f2f6;
        border-radius: 0.5rem;
        padding: 1rem;
        height: 400px;
        overflow-y: auto;
        font-family: monospace;
        margin-bottom: 1rem;
    }

    .footer {
        font-size: 0.8rem;
        color: #888888;
        text-align: center;
        margin-top: 2rem;
    }

    /* Custom CSS for sidebar */
    .css-1d391kg {
        padding-top: 2rem;
    }
    </style>
    """, unsafe_allow_html=True)


def set_page_config():
    """
    Set up the Streamlit page configuration with theme and layout settings.
    """
    st.set_page_config(
        page_title="optiMHC",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )


def main_header(text):
    """
    Display a main header with custom styling.

    Args:
        text: Header text to display
    """
    st.markdown(f'<h1 class="main-header">{text}</h1>', unsafe_allow_html=True)


def sub_header(text):
    """
    Display a sub-header with custom styling.

    Args:
        text: Sub-header text to display
    """
    st.markdown(f'<h2 class="sub-header">{text}</h2>', unsafe_allow_html=True)


def info_box(text):
    """
    Display an information box with custom styling.

    Args:
        text: Information text to display
    """
    st.markdown(f'<div class="info-box">{text}</div>', unsafe_allow_html=True)


def success_box(text):
    """
    Display a success box with custom styling.

    Args:
        text: Success message to display
    """
    st.markdown(f'<div class="success-box">{text}</div>', unsafe_allow_html=True)


def warning_box(text):
    """
    Display a warning box with custom styling.

    Args:
        text: Warning message to display
    """
    st.markdown(f'<div class="warning-box">{text}</div>', unsafe_allow_html=True)


def error_box(text):
    """
    Display an error box with custom styling.

    Args:
        text: Error message to display
    """
    st.markdown(f'<div class="error-box">{text}</div>', unsafe_allow_html=True)


def footer():
    """
    Display a footer with copyright and version information.
    """
    from optimhc import __version__ as version

    st.markdown(
        f'<div class="footer">optiMHC v{version} | '
        f'© {2023}-{2024} Zixiang Shang</div>',
        unsafe_allow_html=True
    )
