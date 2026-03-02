"""
File upload component for optiMHC GUI.
"""

import os
import streamlit as st
import yaml
import tempfile
from typing import Dict, Any, Optional, Tuple
from optimhc.gui.utils import load_config_from_yaml


def config_file_uploader() -> Optional[Dict[str, Any]]:
    """
    Display a file uploader for configuration files.

    Returns:
        Configuration dictionary if a file is uploaded, None otherwise
    """
    uploaded_file = st.file_uploader(
        "Upload configuration file",
        type=["yaml", "yml"],
        help="Upload a YAML configuration file"
    )

    if uploaded_file is not None:
        try:
            config = yaml.safe_load(uploaded_file)
            st.success(f"Configuration file '{uploaded_file.name}' loaded successfully")
            return config
        except Exception as e:
            st.error(f"Error loading configuration file: {str(e)}")
            return None

    return None


def input_path_field(input_type: str, value: str = "", placeholder: str = "") -> str:
    """
    Display an input field for file paths.

    Args:
        input_type: Type of input (pepxml, pin, mzML directory)
        value: Current value
        placeholder: Placeholder text

    Returns:
        String containing file paths, one per line
    """
    if input_type.lower() in ["pepxml", "pin"]:
        help_text = f"Enter the full path to your {input_type} files, one per line"
    elif input_type.lower() == "mzml":
        help_text = "Enter the full path to your mzML directory"
    else:
        help_text = "Enter file paths, one per line"

    paths = st.text_area(
        f"{input_type} File Paths",
        value=value,
        placeholder=placeholder,
        height=100,
        help=help_text,
        key=f"{input_type.lower()}_paths"
    )

    return paths


def yaml_example(example_type: str = "class_i") -> str:
    """
    Return an example YAML configuration by reading from example files.

    Args:
        example_type: Type of example (class_i or class_ii)

    Returns:
        Example YAML configuration as a string
    """
    # Get the path to the examples directory
    current_dir = os.path.dirname(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    )
    examples_dir = os.path.join(current_dir, "examples")

    # Determine which example file to read
    if example_type.lower() == "class_i":
        example_file = os.path.join(examples_dir, "classI_example.yaml")
        default_config = """
experimentName: class_I_example
inputType: pepxml
inputFile:
  - ./data/YE_20180428_SK_HLA_A0202_3Ips_a50mio_R1_01.pep.xml
decoyPrefix: DECOY_
outputDir: ./examples/results/class_I_example
visualization: True
removePreNxtAA: False
numProcesses: 4
showProgress: True
modificationMap:
  "147.035385": "UNIMOD:35"  # Oxidation (M)
  "160.030649": "UNIMOD:4"   # Carbamidomethyl (C)
allele:
  - HLA-A*02:02
featureGenerator:
  - name: Basic
  - name: PWM
    params:
      class: I
  - name: MHCflurry
  - name: NetMHCpan
rescore:
  testFDR: 0.01
  model: Percolator
  numJobs: 4
"""
    else:  # class_ii
        example_file = os.path.join(examples_dir, "classII_example.yaml")
        default_config = """
experimentName: class_II_example
inputType: pepxml
inputFile:
  - ./data/AG20201214_FAIMS_DPB0101_DPA0201_93e6_1hr.pep.xml
decoyPrefix: DECOY_
outputDir: ./examples/results/class_II_example
visualization: True
removePreNxtAA: False
numProcesses: 4
showProgress: True
modificationMap:
  "147.035385": "UNIMOD:35"  # Oxidation (M)
  "160.030649": "UNIMOD:4"   # Carbamidomethyl (C)
allele:
  - HLA-DPA1*02:01-DPB1*01:01
featureGenerator:
  - name: Basic
  - name: PWM
    params:
      class: II
  - name: NetMHCIIpan
rescore:
  testFDR: 0.01
  model: Percolator
  numJobs: 4
"""

    try:
        if not os.path.exists(example_file):
            st.warning(f"Example file not found: {example_file}, using default configuration")
            return default_config

        with open(example_file, 'r') as f:
            content = f.read()
            if not content:
                st.warning(f"Example file is empty: {example_file}, using default configuration")
                return default_config
            return content
    except Exception as e:
        st.warning(f"Error reading example file: {str(e)}, using default configuration")
        return default_config
