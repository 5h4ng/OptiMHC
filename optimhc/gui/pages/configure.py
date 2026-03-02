"""
Configuration page for the optiMHC GUI.
"""

import streamlit as st
import os
import yaml
from typing import Dict, Any

from optimhc.gui.style import main_header, sub_header, info_box, warning_box
from optimhc.gui.components.config_form import validate_config, render_config_summary
from optimhc.gui.components.file_upload import config_file_uploader, input_path_field, yaml_example
from optimhc.gui.components.feature_generator_form import feature_generator_form
from optimhc.gui.components.rescore_form import rescore_form
from optimhc.gui.components.modification_map import modification_map_form


def render():
    """
    Render the configuration page.
    """
    main_header("Configure Pipeline")

    st.markdown("""
    Configure your optiMHC pipeline by filling out the form below, or upload an existing configuration file.
    """)

    # Initialize session state for configuration
    if "config" not in st.session_state:
        st.session_state.config = {}

    # File upload option and example buttons
    col1, col2, col3 = st.columns([2, 1, 1])

    with col1:
        uploaded_config = config_file_uploader()

        if uploaded_config:
            # Validate configuration
            errors = validate_config(uploaded_config)

            if errors:
                for error in errors:
                    st.error(error)
            else:
                # Save configuration to session state
                st.session_state.config = uploaded_config

                # Handle modification map in session state for the form
                if "modificationMap" in uploaded_config and uploaded_config["modificationMap"]:
                    st.session_state.num_modifications = len(uploaded_config["modificationMap"])
                    st.session_state.modification_masses = list(uploaded_config["modificationMap"].keys())
                    st.session_state.modification_values = list(uploaded_config["modificationMap"].values())

                st.success("Configuration loaded successfully!")
                st.rerun()

    with col2:
        if st.button("Load Class I Example", use_container_width=True):
            example_config = yaml.safe_load(yaml_example("class_i"))
            # Ensure modification map is included
            if "modificationMap" not in example_config or not example_config["modificationMap"]:
                example_config["modificationMap"] = {
                    "147.035385": "UNIMOD:35",  # Oxidation (M) - full modified residue mass
                    "160.030649": "UNIMOD:4"    # Carbamidomethyl (C) - full modified residue mass
                }
            st.session_state.config = example_config
            # Also update the modification map in session state for the form
            st.session_state.num_modifications = len(example_config["modificationMap"])
            st.session_state.modification_masses = list(example_config["modificationMap"].keys())
            st.session_state.modification_values = list(example_config["modificationMap"].values())
            st.success("Class I example configuration loaded!")
            st.rerun()

    with col3:
        if st.button("Load Class II Example", use_container_width=True):
            example_config = yaml.safe_load(yaml_example("class_ii"))
            # Ensure modification map is included
            if "modificationMap" not in example_config or not example_config["modificationMap"]:
                example_config["modificationMap"] = {
                    "147.035385": "UNIMOD:35",  # Oxidation (M) - full modified residue mass
                    "160.030649": "UNIMOD:4"    # Carbamidomethyl (C) - full modified residue mass
                }
            st.session_state.config = example_config
            # Also update the modification map in session state for the form
            st.session_state.num_modifications = len(example_config["modificationMap"])
            st.session_state.modification_masses = list(example_config["modificationMap"].keys())
            st.session_state.modification_values = list(example_config["modificationMap"].values())
            st.success("Class II example configuration loaded!")
            st.rerun()

    # Setup configuration tabs
    tabs = st.tabs(["Basic Settings", "Feature Generators", "Rescoring", "Summary"])

    # Tab 1: Basic Settings
    with tabs[0]:
        st.subheader("Basic Settings")

        # Experiment name
        experiment_name = st.text_input(
            "Experiment Name",
            value=st.session_state.config.get("experimentName", ""),
            placeholder="my_experiment",
            help="Name of the experiment"
        )

        # Input type and files
        col1, col2 = st.columns(2)
        with col1:
            input_type = st.selectbox(
                "Input Type",
                options=["pepxml", "pin"],
                index=["pepxml", "pin"].index(st.session_state.config.get("inputType", "pepxml")),
                help="Type of input file"
            )

        with col2:
            decoy_prefix = st.text_input(
                "Decoy Prefix",
                value=st.session_state.config.get("decoyPrefix", "DECOY_"),
                help="Prefix used to identify decoy sequences"
            )

        # Input files - use direct path input
        input_files = st.session_state.config.get("inputFile", [])
        if isinstance(input_files, str):
            input_files = [input_files]

        input_files_str = input_path_field(
            input_type,
            value="\n".join(input_files) if input_files else "",
            placeholder="../examples/data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml"
        )

        # Output directory
        output_dir = st.text_input(
            "Output Directory",
            value=st.session_state.config.get("outputDir", "../examples/results"),
            placeholder="../examples/results/my_experiment",
            help="Directory where results will be saved"
        )

        # Allele settings
        st.subheader("Allele Settings")

        alleles = st.session_state.config.get("allele", [])
        if isinstance(alleles, str):
            alleles = [alleles]

        alleles_str = st.text_area(
            "Alleles",
            value="\n".join(alleles) if alleles else "",
            placeholder="HLA-A*02:01\nHLA-B*07:02",
            height=100,
            help="One allele per line, e.g., HLA-A*02:01 for class I or HLA-DPA1*02:01-DPB1*01:01 for class II"
        )

        # Modification map
        st.subheader("Modification Map")
        st.markdown("Maps full modified residue masses to UNIMOD identifiers.")
        modification_map = modification_map_form(st.session_state.config.get("modificationMap", {}))

        # Performance settings
        st.subheader("Performance Settings")

        col1, col2 = st.columns(2)

        with col1:
            num_processes = st.number_input(
                "Number of Processes",
                min_value=1,
                max_value=64,
                value=int(st.session_state.config.get("numProcesses", 4)),
                help="Number of parallel processes"
            )

        with col2:
            show_progress = st.checkbox(
                "Show Progress",
                value=st.session_state.config.get("showProgress", True),
                help="Show progress bars during processing"
            )

        col1, col2 = st.columns(2)

        with col1:
            visualization = st.checkbox(
                "Enable Visualization",
                value=st.session_state.config.get("visualization", True),
                help="Generate visualizations of results"
            )

        with col2:
            remove_pre_nxt_aa = st.checkbox(
                "Remove Pre/Next Amino Acids",
                value=st.session_state.config.get("removePreNxtAA", False),
                help="Remove pre/post neighboring amino acids in sequence processing"
            )

        log_level = st.selectbox(
            "Log Level",
            options=["DEBUG", "INFO", "WARNING", "ERROR"],
            index=["DEBUG", "INFO", "WARNING", "ERROR"].index(st.session_state.config.get("logLevel", "INFO")),
            help="Logging verbosity level"
        )

        # Save basic settings
        if st.button("Save Basic Settings", key="save_basic"):
            # Process input fields
            input_files = [f for f in input_files_str.strip().split("\n") if f]
            alleles = [a for a in alleles_str.strip().split("\n") if a]

            # Update configuration
            st.session_state.config.update({
                "experimentName": experiment_name,
                "inputType": input_type,
                "inputFile": input_files,
                "decoyPrefix": decoy_prefix,
                "outputDir": output_dir,
                "allele": alleles,
                "modificationMap": modification_map,
                "numProcesses": num_processes,
                "showProgress": show_progress,
                "visualization": visualization,
                "removePreNxtAA": remove_pre_nxt_aa,
                "logLevel": log_level
            })

            st.success("Basic settings saved!")
            st.rerun()

    # Tab 2: Feature Generators
    with tabs[1]:
        # Feature generators configuration
        feature_generators = feature_generator_form(st.session_state.config.get("featureGenerator", []))

        # Save feature generators
        if st.button("Save Feature Generators", key="save_features"):
            st.session_state.config["featureGenerator"] = feature_generators
            st.success("Feature generators saved!")

    # Tab 3: Rescoring
    with tabs[2]:
        # Rescoring configuration
        rescore_config = rescore_form(st.session_state.config.get("rescore", {}))

        # Save rescoring settings
        if st.button("Save Rescoring Settings", key="save_rescore"):
            st.session_state.config["rescore"] = rescore_config
            st.success("Rescoring settings saved!")

    # Tab 4: Summary
    with tabs[3]:
        # Full configuration summary
        if st.session_state.config:
            st.subheader("Configuration Summary")

            # Display configuration summary
            render_config_summary(st.session_state.config)

            # Validate configuration
            errors = validate_config(st.session_state.config)

            if errors:
                st.error("Configuration has errors:")
                for error in errors:
                    st.error(f"- {error}")
            else:
                st.success("Configuration is valid!")

            # Option to download configuration
            if st.button("Download Configuration as YAML"):
                # Create YAML content
                config_yaml = yaml.dump(st.session_state.config, default_flow_style=False)

                # Use streamlit's download button
                filename = f"{st.session_state.config.get('experimentName', 'optimhc_config')}.yaml"
                st.download_button(
                    label="Download Configuration File",
                    data=config_yaml,
                    file_name=filename,
                    mime="text/yaml"
                )
        else:
            st.warning("No configuration has been saved yet. Please fill out and save the forms in the other tabs.")

    # Navigation buttons
    st.markdown("---")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("← Back to Home"):
            st.session_state.page = "home"
            st.rerun()

    with col2:
        if st.button("Continue to Run →"):
            # Check if we have a valid configuration
            if not st.session_state.config:
                warning_box("Please configure the pipeline first.")
            else:
                errors = validate_config(st.session_state.config)
                if errors:
                    for error in errors:
                        st.error(error)
                else:
                    st.session_state.page = "run"
                    st.rerun()
