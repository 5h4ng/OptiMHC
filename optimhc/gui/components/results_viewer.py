"""
Results viewer component for optiMHC GUI.
"""

import os
import glob
import base64
from typing import List, Dict, Any, Optional
import streamlit as st
import pandas as pd
import plotly.express as px
from pathlib import Path

from optimhc.gui.utils import scan_output_directory


def get_image_as_base64(file_path):
    """
    Get image file as base64 string.

    Args:
        file_path: Path to image file

    Returns:
        Base64 encoded image
    """
    with open(file_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode()


def display_image(file_path, caption=None):
    """
    Display an image with caption.

    Args:
        file_path: Path to image file
        caption: Optional caption for the image
    """
    try:
        # Use HTML to have more control over image sizing
        img_format = file_path.split('.')[-1].lower()
        img_base64 = get_image_as_base64(file_path)
        html = f'<img src="data:image/{img_format};base64,{img_base64}" style="max-width: 100%; height: auto;">'

        if caption:
            html = f'{html}<div style="text-align: center; font-style: italic; margin-top: 5px;">{caption}</div>'

        st.markdown(html, unsafe_allow_html=True)
    except Exception as e:
        st.error(f"Error displaying image {os.path.basename(file_path)}: {str(e)}")


def display_csv(file_path, caption=None):
    """
    Display a CSV file as a dataframe.

    Args:
        file_path: Path to CSV file
        caption: Optional caption for the table
    """
    try:
        # Determine file type and use appropriate reader
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path)
        elif file_path.endswith('.tsv'):
            df = pd.read_csv(file_path, sep='\t')
        elif file_path.endswith('.xlsx'):
            df = pd.read_excel(file_path)
        else:
            st.warning(f"Unsupported file format: {os.path.basename(file_path)}")
            return

        if caption:
            st.markdown(f"**{caption}**")

        # Display dataframe with pagination if large
        if len(df) > 100:
            # Add pagination controls
            page_size = st.slider(f"Rows per page for {os.path.basename(file_path)}", 10, 100, 50)
            page = st.number_input(f"Page for {os.path.basename(file_path)}", min_value=1, max_value=max(1, len(df) // page_size + 1), value=1)

            start_idx = (page - 1) * page_size
            end_idx = min(start_idx + page_size, len(df))

            st.dataframe(df.iloc[start_idx:end_idx])
            st.text(f"Showing rows {start_idx+1}-{end_idx} of {len(df)}")
        else:
            st.dataframe(df)

        # Option to download the file
        csv_data = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label=f"Download {os.path.basename(file_path)}",
            data=csv_data,
            file_name=os.path.basename(file_path),
            mime="text/csv"
        )
    except Exception as e:
        st.error(f"Error displaying data from {os.path.basename(file_path)}: {str(e)}")


def results_viewer(output_dir: str):
    """
    Display results from the output directory, focusing only on figures.

    Args:
        output_dir: Path to the output directory
    """
    st.subheader("Results Visualization")

    if not os.path.exists(output_dir):
        st.warning(f"Output directory '{output_dir}' does not exist yet")
        return

    # Check if there's an experiment name in the session state
    experiment_name = None
    if "config" in st.session_state:
        experiment_name = st.session_state.config.get("experimentName")

    # Look for figures in experiment subdirectory if it exists
    if experiment_name:
        experiment_dir = os.path.join(output_dir, experiment_name)
        if os.path.exists(experiment_dir):
            # Update output_dir to use the experiment directory
            output_dir = experiment_dir

    # Check for figures subfolder
    figures_dir = os.path.join(output_dir, "figures")
    has_figures_subdir = os.path.exists(figures_dir) and os.path.isdir(figures_dir)

    # Get all figure files
    figure_files = []

    # First check figures subdirectory if it exists
    if has_figures_subdir:
        figure_files = glob.glob(os.path.join(figures_dir, "*.png"))
        figure_files.extend(glob.glob(os.path.join(figures_dir, "*.jpg")))
        figure_files.extend(glob.glob(os.path.join(figures_dir, "*.svg")))
        figure_files.extend(glob.glob(os.path.join(figures_dir, "*.pdf")))

    # If no figures found in subfolder (or none exists), check main directory
    if not figure_files:
        figure_files = glob.glob(os.path.join(output_dir, "*.png"))
        figure_files.extend(glob.glob(os.path.join(output_dir, "*.jpg")))
        figure_files.extend(glob.glob(os.path.join(output_dir, "*.svg")))
        figure_files.extend(glob.glob(os.path.join(output_dir, "*.pdf")))

    # Show figures
    if figure_files:
        figure_count = len(figure_files)

        # Show information about found figures
        st.caption(f"Found {figure_count} figures in {figures_dir if has_figures_subdir else output_dir}")

        # Create a selector for figures if there are more than one
        if figure_count > 1:
            # Sort figures by name for consistency
            figure_files.sort()
            figure_names = [os.path.basename(f) for f in figure_files]

            # Figure selection dropdown
            selected_figure = st.selectbox(
                "Select figure to view:",
                figure_names,
                key="figure_selector"
            )

            selected_index = figure_names.index(selected_figure)

            # Display the selected figure
            st.markdown("### Selected Figure")
            display_image(figure_files[selected_index], caption=selected_figure)

            # Show thumbnail gallery with 3 columns
            st.markdown("### All Figures")

            # Use columns for gallery display
            cols = st.columns(min(3, figure_count))
            for i, file_path in enumerate(figure_files):
                with cols[i % 3]:
                    try:
                        img_format = file_path.split('.')[-1].lower()
                        if img_format in ['png', 'jpg', 'jpeg', 'svg']:
                            img_base64 = get_image_as_base64(file_path)
                            html = f'<img src="data:image/{img_format};base64,{img_base64}" style="width: 100%; height: auto;">'
                            st.markdown(html, unsafe_allow_html=True)
                            st.caption(os.path.basename(file_path))
                        else:
                            # For non-image formats like PDF
                            st.info(f"{os.path.basename(file_path)} (PDF file)")
                    except Exception as e:
                        st.error(f"Error with image {os.path.basename(file_path)}")
        else:
            # Just show the single figure
            display_image(figure_files[0], caption=os.path.basename(figure_files[0]))
    else:
        st.info(f"No figures found in the output directory. Run the pipeline to generate results.")

    # Refresh button
    if st.button("Refresh Results", key="refresh_results_btn"):
        st.rerun()


def results_summary(output_dir: str):
    """
    Display a summary of results.

    Args:
        output_dir: Path to the output directory
    """
    if not os.path.exists(output_dir):
        return

    # Look for summary statistics
    summary_files = glob.glob(os.path.join(output_dir, "**/summary*.csv"), recursive=True)
    summary_files.extend(glob.glob(os.path.join(output_dir, "**/stats*.csv"), recursive=True))

    if not summary_files:
        return

    st.subheader("Results Summary")

    try:
        # Use the first summary file found
        df = pd.read_csv(summary_files[0])

        # Display key statistics
        cols = st.columns(3)

        if 'total_psms' in df.columns:
            with cols[0]:
                st.metric("Total PSMs", df['total_psms'].iloc[0])

        if 'target_psms' in df.columns:
            with cols[1]:
                st.metric("Target PSMs", df['target_psms'].iloc[0])

        if 'decoy_psms' in df.columns:
            with cols[2]:
                st.metric("Decoy PSMs", df['decoy_psms'].iloc[0])

        # Generate a simple chart if possible
        if {'fdr', 'target_psms'}.issubset(df.columns):
            st.subheader("FDR vs. Target PSMs")
            fig = px.line(df, x='fdr', y='target_psms', title='PSMs at different FDR thresholds')
            st.plotly_chart(fig, use_container_width=True)
    except Exception as e:
        st.error(f"Error generating results summary: {str(e)}")
