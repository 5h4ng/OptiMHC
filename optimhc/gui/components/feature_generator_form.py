"""
Feature generator form component for optiMHC GUI.
"""

import streamlit as st
from typing import Dict, Any, List

def feature_generator_form(existing_generators: List[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """
    Create a form for configuring feature generators.

    Args:
        existing_generators: List of existing feature generator configurations

    Returns:
        List of feature generator configurations
    """
    if existing_generators is None:
        existing_generators = []

    # Convert existing generators to a dict for easier lookup
    existing_gen_dict = {}
    for gen in existing_generators:
        existing_gen_dict[gen["name"]] = gen.get("params", {})

    feature_generators = []

    st.subheader("Feature Generators")

    # Determine the MHC class from existing configuration
    # Look for PWM first as it has explicit class parameter
    mhc_class = None
    for gen in existing_generators:
        if gen["name"] == "PWM" and "params" in gen and "class" in gen["params"]:
            mhc_class = gen["params"]["class"]
            break

    # If PWM not found, infer from presence of NetMHCIIpan
    if mhc_class is None:
        if any(gen["name"] == "NetMHCIIpan" for gen in existing_generators):
            mhc_class = "II"
        else:
            mhc_class = "I"  # Default to class I

    # MHC class selection
    mhc_class = st.radio(
        "MHC Class",
        options=["I", "II"],
        index=0 if mhc_class == "I" else 1,
        horizontal=True,
        help="Select MHC class for appropriate feature generators."
    )

    st.markdown("---")
    st.markdown("Select which feature generators to use in the pipeline:")

    # Basic feature generator (always available)
    if st.checkbox("Basic", value="Basic" in existing_gen_dict or not existing_generators, key="basic_gen"):
        feature_generators.append({"name": "Basic"})

    # PWM feature generator (class parameter set automatically based on MHC class selection)
    if st.checkbox("PWM", value="PWM" in existing_gen_dict, key="pwm_gen"):
        feature_generators.append({
            "name": "PWM",
            "params": {"class": mhc_class}
        })

    # Class I specific generators
    if mhc_class == "I":
        # MHCflurry (class I only)
        if st.checkbox("MHCflurry", value="MHCflurry" in existing_gen_dict, key="mhcflurry_gen"):
            feature_generators.append({"name": "MHCflurry"})

        # NetMHCpan (class I only)
        if st.checkbox("NetMHCpan", value="NetMHCpan" in existing_gen_dict, key="netmhcpan_gen"):
            feature_generators.append({"name": "NetMHCpan"})

    # Class II specific generators
    else:  # mhc_class == "II"
        # NetMHCIIpan (class II only)
        if st.checkbox("NetMHCIIpan", value="NetMHCIIpan" in existing_gen_dict, key="netmhciipan_gen"):
            feature_generators.append({"name": "NetMHCIIpan"})

    # DeepLC feature generator (available for both classes)
    if st.checkbox("DeepLC", value="DeepLC" in existing_gen_dict, key="deeplc_gen"):
        deeplc_params = {}

        col1, col2 = st.columns(2)
        with col1:
            calibration_criteria = st.text_input(
                "Calibration Criteria",
                value=existing_gen_dict.get("DeepLC", {}).get("calibrationCriteria", "expect"),
                key="deeplc_calibration_criteria",
                help="Criteria for calibration (e.g., expect, xcorr, hyperscore)"
            )
            deeplc_params["calibrationCriteria"] = calibration_criteria

        with col2:
            lower_is_better = st.checkbox(
                "Lower Is Better",
                value=existing_gen_dict.get("DeepLC", {}).get("lowerIsBetter", True),
                key="deeplc_lower_is_better",
                help="Whether lower values of the calibration criteria are better (True for expect, False for xcorr/hyperscore)"
            )
            deeplc_params["lowerIsBetter"] = lower_is_better

        calibration_size = st.slider(
            "Calibration Size",
            min_value=0.05,
            max_value=0.5,
            value=float(existing_gen_dict.get("DeepLC", {}).get("calibrationSize", 0.1)),
            step=0.05,
            key="deeplc_calibration_size",
            help="Fraction of PSMs to use for calibration (0.05-0.5)"
        )
        deeplc_params["calibrationSize"] = calibration_size

        feature_generators.append({"name": "DeepLC", "params": deeplc_params})

    # SpectralSimilarity feature generator (with AlphaPeptDeep as default)
    if st.checkbox("SpectralSimilarity", value="SpectralSimilarity" in existing_gen_dict or not existing_generators, key="spectra_similarity_gen"):
        ss_params = {}

        st.markdown("#### SpectralSimilarity Settings")

        model = st.selectbox(
            "Model",
            options=["AlphaPeptDeep_ms2_generic"],
            index=["AlphaPeptDeep_ms2_generic"].index(
                existing_gen_dict.get("SpectralSimilarity", {}).get("model", "AlphaPeptDeep_ms2_generic")
            ),
            key="spectra_similarity_model",
            help="Prediction model for theoretical spectra"
        )
        ss_params["model"] = model

        instrument = st.selectbox(
            "Instrument",
            options=["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"],
            index=["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"].index(
                existing_gen_dict.get("SpectralSimilarity", {}).get("instrument", "LUMOS")
            ),
            key="spectra_similarity_instrument",
            help="Available instruments: QE, LUMOS, TIMSTOF, SCIEXTOF"
        )
        ss_params["instrument"] = instrument

        # mzML directory path
        mzml_dir = st.text_input(
            "mzML Directory Path",
            value=existing_gen_dict.get("SpectralSimilarity", {}).get("mzmlDir", "./data"),
            key="spectra_similarity_mzml_dir",
            help="Path to directory containing mzML files"
        )
        if mzml_dir:
            ss_params["mzmlDir"] = mzml_dir

        # Spectrum ID pattern
        spectrum_id_pattern = st.text_input(
            "Spectrum ID Pattern",
            value=existing_gen_dict.get("SpectralSimilarity", {}).get("spectrumIdPattern", "(.+?)\\.\\d+\\.\\d+\\.\\d+"),
            key="spectra_similarity_spectrum_id_pattern",
            help="Regular expression pattern to extract mzML filename from spectrum IDs. Default pattern: (.+?)\\.\\d+\\.\\d+\\.\\d+"
        )
        if spectrum_id_pattern:
            ss_params["spectrumIdPattern"] = spectrum_id_pattern

        collision_energy = st.number_input(
            "Collision Energy",
            min_value=20,
            max_value=40,
            value=int(existing_gen_dict.get("SpectralSimilarity", {}).get("collisionEnergy", 28)),
            key="spectra_similarity_collision_energy",
            help="Collision energy used during acquisition (typical range: 25-30)"
        )
        ss_params["collisionEnergy"] = collision_energy

        tolerance = st.slider(
            "Tolerance (ppm)",
            min_value=10,
            max_value=50,
            value=int(existing_gen_dict.get("SpectralSimilarity", {}).get("tolerance", 20)),
            step=5,
            key="spectra_similarity_tolerance",
            help="Mass tolerance in ppm for peak matching (10-50 ppm)"
        )
        ss_params["tolerance"] = tolerance

        num_top_peaks = st.slider(
            "Number of Top Peaks",
            min_value=10,
            max_value=100,
            value=int(existing_gen_dict.get("SpectralSimilarity", {}).get("numTopPeaks", 36)),
            step=2,
            key="spectra_similarity_num_top_peaks",
            help="Number of most intense peaks to consider for matching"
        )
        ss_params["numTopPeaks"] = num_top_peaks

        url = st.text_input(
            "API URL",
            value=existing_gen_dict.get("SpectralSimilarity", {}).get("url", "koina.wilhelmlab.org:443"),
            key="spectra_similarity_url",
            help="AlphaPept API URL (default: koina.wilhelmlab.org:443)"
        )
        if url:
            ss_params["url"] = url

        feature_generators.append({"name": "SpectralSimilarity", "params": ss_params})

    # OverlappingPeptide feature generator
    if st.checkbox("OverlappingPeptide", value="OverlappingPeptide" in existing_gen_dict, key="overlapping_peptide_gen"):
        op_params = {}

        st.markdown("#### OverlappingPeptide Settings")

        col1, col2 = st.columns(2)
        with col1:
            min_overlap_length = st.number_input(
                "Min Overlap Length",
                min_value=5,
                max_value=15,
                value=int(existing_gen_dict.get("OverlappingPeptide", {}).get("minOverlapLength", 7)),
                key="op_min_overlap_length",
                help="Minimum number of amino acids that must overlap"
            )
            op_params["minOverlapLength"] = min_overlap_length

        with col2:
            overlapping_score = st.text_input(
                "Overlapping Score",
                value=existing_gen_dict.get("OverlappingPeptide", {}).get("overlappingScore", "expect"),
                key="op_overlapping_score",
                help="Score to use for overlapping peptides (e.g., expect, xcorr, hyperscore)"
            )
            op_params["overlappingScore"] = overlapping_score

        col1, col2 = st.columns(2)
        with col1:
            min_length = st.number_input(
                "Min Length",
                min_value=5,
                max_value=15,
                value=int(existing_gen_dict.get("OverlappingPeptide", {}).get("minLength", 7 if mhc_class == "I" else 9)),
                key="op_min_length",
                help="Minimum peptide length to consider"
            )
            op_params["minLength"] = min_length

        with col2:
            max_length = st.number_input(
                "Max Length",
                min_value=10,
                max_value=50,
                value=int(existing_gen_dict.get("OverlappingPeptide", {}).get("maxLength", 20 if mhc_class == "I" else 30)),
                key="op_max_length",
                help="Maximum peptide length to consider"
            )
            op_params["maxLength"] = max_length

        feature_generators.append({"name": "OverlappingPeptide", "params": op_params})

    # Warning if no generators selected
    if not feature_generators:
        st.warning("Please select at least one feature generator.")

    return feature_generators
