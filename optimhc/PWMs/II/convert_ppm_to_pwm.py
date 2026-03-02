#!/usr/bin/env python3

import math
import os

# Define the background frequency
BACKGROUND_FREQ = 1.0 / 20  # 0.05


def ppm_to_pwm(ppm_matrix):
    """
    Convert a PPM matrix to a PWM matrix using log2(PPM / background_freq).
    """
    pwm_matrix = {}
    for aa, freqs in ppm_matrix.items():
        pwm_values = []
        for p in freqs:
            if p <= 0:
                raise ValueError(
                    f"PPM value for amino acid '{aa}' is zero or negative, cannot compute log."
                )
            pwm = math.log2(p / BACKGROUND_FREQ)
            pwm_values.append(pwm)
        pwm_matrix[aa] = pwm_values
    return pwm_matrix


def read_ppm_file(filepath):
    """
    Read a PPM file and return a dictionary with amino acids as keys and lists of frequencies as values.
    """
    ppm = {}
    with open(filepath, "r") as f:
        for line in f:
            if line.strip() == "":
                continue  # Skip empty lines
            parts = line.strip().split()
            aa = parts[0]
            freqs = [float(x) for x in parts[1:]]
            ppm[aa] = freqs
    return ppm


def write_pwm_file(filepath, pwm_matrix):
    """
    Write the PWM matrix to a file, overwriting the original PPM file.
    """
    with open(filepath, "w") as f:
        for aa in sorted(pwm_matrix.keys()):
            pwm_values = pwm_matrix[aa]
            # Format each PWM value to six decimal places
            formatted_values = " ".join(f"{v:.6f}" for v in pwm_values)
            f.write(f"{aa}\t{formatted_values}\n")


def process_file(filepath):
    """
    Process a single PPM file: convert it to PWM and overwrite the original file.
    """
    try:
        ppm = read_ppm_file(filepath)
        pwm = ppm_to_pwm(ppm)
        write_pwm_file(filepath, pwm)
        print(f"Converted: {filepath}")
    except Exception as e:
        print(f"Error processing {filepath}: {e}")


def main():
    """
    Main function to traverse the II directory and convert all PPM files to PWM.
    """
    base_dir = os.getcwd()

    # Process c_flank_ppm and n_flank_ppm in the base directory
    for filename in ["c_flank_ppm", "n_flank_ppm"]:
        filepath = os.path.join(base_dir, filename)
        if os.path.isfile(filepath):
            process_file(filepath)
        else:
            print(f"File not found: {filepath}")

    # Traverse all subdirectories (alleles) and process their PPM files
    for entry in os.listdir(base_dir):
        subdir = os.path.join(base_dir, entry)
        if os.path.isdir(subdir):
            # Assuming there's only one PPM file per allele directory with a .txt extension
            for subfile in os.listdir(subdir):
                if subfile.endswith(".txt"):
                    subfilepath = os.path.join(subdir, subfile)
                    process_file(subfilepath)


if __name__ == "__main__":
    main()
