"""
Modification map component for optiMHC GUI.
"""

import streamlit as st
from typing import Dict, Any, Optional

def modification_map_form(existing_map: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """
    Create a form for modification map configuration.

    Args:
        existing_map: Existing modification map configuration

    Returns:
        Modification map dictionary mapping masses to UNIMOD values
    """
    if existing_map is None:
        existing_map = {
            "147.035385": "UNIMOD:35",  # Oxidation (M) - Full modified residue mass
            "160.030649": "UNIMOD:4",   # Carbamidomethyl (C) - Full modified residue mass
            "166.998359": "UNIMOD:21"   # Phospho (S) - Full modified residue mass
        }

    st.subheader("Modification Map")

    st.markdown("""
    Specify the mapping from modification masses to UNIMOD identifiers.
    The mass value should be the FULL modified residue mass (amino acid + modification) as found in pepXML parameters.
    All modifications need to be explicitly encoded in the sequence (e.g., C[UNIMOD:4] for carbamidomethylated cysteine).
    """)

    # Create a container for the dynamic map
    modification_map = {}

    # Use session state to track number of modification entries
    if "num_modifications" not in st.session_state:
        st.session_state.num_modifications = len(existing_map)
        st.session_state.modification_masses = list(existing_map.keys())
        st.session_state.modification_values = list(existing_map.values())

    # Add/remove modification controls
    col1, col2 = st.columns([1, 5])
    with col1:
        if st.button("➕ Add Modification", key="add_modification"):
            st.session_state.num_modifications += 1
            st.session_state.modification_masses.append("")
            st.session_state.modification_values.append("UNIMOD:")
            st.rerun()
    with col2:
        if st.session_state.num_modifications > 0 and st.button("➖ Remove Last Modification", key="remove_modification"):
            st.session_state.num_modifications -= 1
            if st.session_state.modification_masses:
                st.session_state.modification_masses.pop()
            if st.session_state.modification_values:
                st.session_state.modification_values.pop()
            st.rerun()

    # Create a table-like interface for modifications
    if st.session_state.num_modifications > 0:
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Mass (Residue+Modification)**")
        with col2:
            st.markdown("**UNIMOD Identifier**")

        for i in range(st.session_state.num_modifications):
            col1, col2 = st.columns(2)
            with col1:
                mass = st.text_input(
                    "Mass",
                    value=st.session_state.modification_masses[i] if i < len(st.session_state.modification_masses) else "",
                    key=f"mod_mass_{i}",
                    label_visibility="collapsed"
                )
                st.session_state.modification_masses[i] = mass

            with col2:
                unimod = st.text_input(
                    "UNIMOD",
                    value=st.session_state.modification_values[i] if i < len(st.session_state.modification_values) else "UNIMOD:",
                    key=f"mod_unimod_{i}",
                    label_visibility="collapsed"
                )
                st.session_state.modification_values[i] = unimod

            # Add to modification map
            if mass and unimod:
                modification_map[mass] = unimod

    # Information about common modifications
    with st.expander("Common Modifications (Note: Values are examples, check your pepXML)", expanded=False):
        st.markdown("""
        | Mass (Full) | UNIMOD ID | Modification | Target Residues |
        |------|-----------|--------------|--------------|
        | 147.035385 | UNIMOD:35 | Oxidation | M |
        | 160.030649 | UNIMOD:4 | Carbamidomethyl | C |

        Note: These are full masses (amino acid + modification). You must check your pepXML file parameters to find the exact masses used in your data.
        """)

    return modification_map
