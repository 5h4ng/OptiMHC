"""
feature_generation.py

Implements feature generation logic for optiMHC, supporting multiple feature generators
(Basic, OverlappingPeptide, PWM, MHCflurry, NetMHCpan, NetMHCIIpan, DeepLC, SpectralSimilarity, etc.).
"""
import os
import logging
import re
import gc

# The reason why we need to import the feature generators here is that 
# the package 'mhctools' affect the logging configuration of optiMHC.
# TODO: find a better way to handle this.
from optimhc.feature_generator.netMHCpan import NetMHCpanFeatureGenerator
from optimhc.feature_generator.netMHCIIpan import NetMHCIIpanFeatureGenerator

logger = logging.getLogger(__name__)


# TODO: refactor the code to pass config as a parameter to the generators
def generate_features(psms, config):
    """
    Generate features from different generators according to the configuration.

    Parameters
    ----------
    psms : PsmContainer
        A container object holding PSMs and relevant data.
    config : dict
        Configuration dictionary loaded from YAML or CLI.

    Returns
    -------
    None
        Features are added in-place to the PsmContainer.

    Examples
    --------
    >>> generate_features(psms, config)
    """
    remove_modification = True
    remove_pre_nxt_aa = config["removePreNxtAA"]
    n_processes = config["numProcesses"]
    show_progress = config["showProgress"]
    output_dir = config["outputDir"]
    mod_dict = config.get("modificationMap", None)
    if mod_dict == {}:
        mod_dict = None
    feature_generators = config.get("featureGenerator", None)
    allele = config.get("allele", None)
    unique_peptides = list(set(psms.peptides))

    if feature_generators is not None:
        for generator_config in feature_generators:
            if not isinstance(generator_config, dict):
                logger.warning(
                    "Feature generator config is not a dictionary, skipping..."
                )
                continue

            generator_type = generator_config.get("name")
            logger.info(f"Generating features with {generator_type}...")
            generator_params = generator_config.get("params", {})

            if generator_type == "OverlappingPeptide":
                from optimhc.feature_generator.overlapping_peptide import (
                    OverlappingPeptideFeatureGenerator,
                    assign_brother_aggregated_feature,
                )
                overlapping_peptide = OverlappingPeptideFeatureGenerator(
                    unique_peptides,
                    min_overlap_length=generator_params.get("minOverlapLength", 8),
                    min_length=generator_params.get("minLength", 8),
                    max_length=generator_params.get("maxLength", 25),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    remove_modification=remove_modification,
                )
                overlapping_features = overlapping_peptide.generate_features()
                full_data = overlapping_peptide.get_full_data()

                psms.add_metadata(
                    full_data[["Peptide", "contig_member_count", "ContigSequence"]],
                    psms_key=psms.peptide_column,
                    metadata_key="Peptide",
                    source="OverlappingPeptide",
                )
                psms.add_features(
                    overlapping_features,
                    psms_key=psms.peptide_column,
                    feature_key=overlapping_peptide.id_column,
                    source="OverlappingPeptide",
                )
                score = generator_params.get("overlappingScore", None)
                if score:
                    assign_brother_aggregated_feature(
                        psms,
                        feature_columns=score,
                        overlapping_source="OverlappingPeptide",
                        source_name="ContigFeatures",
                    )

                del overlapping_peptide, overlapping_features, full_data
                gc.collect()

            elif generator_type == "Basic":
                from optimhc.feature_generator.basic import BasicFeatureGenerator
                basic_generator = BasicFeatureGenerator(
                    psms.psms[psms.peptide_column].tolist(),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    remove_modification=remove_modification,
                )
                basic_features = basic_generator.generate_features()
                psms.add_features_by_index(
                    basic_features[basic_generator.feature_columns], source="Basic"
                )

                del basic_generator, basic_features
                gc.collect()

            elif generator_type == "PWM":
                from optimhc.feature_generator.PWM import PWMFeatureGenerator
                pwm_generator = PWMFeatureGenerator(
                    unique_peptides,
                    alleles=allele,
                    mhc_class=generator_params.get("class", "I"),
                    remove_modification=remove_modification,
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                )
                pwm_features = pwm_generator.generate_features()
                psms.add_features(
                    pwm_features,
                    psms_key=psms.peptide_column,
                    feature_key=pwm_generator.id_column,
                    source="PWM",
                )

                del pwm_generator, pwm_features
                gc.collect()

            elif generator_type == "MHCflurry":
                from optimhc.feature_generator.mhcflurry import MHCflurryFeatureGenerator
                mhcflurry_generator = MHCflurryFeatureGenerator(
                    unique_peptides,
                    alleles=allele,
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    remove_modification=remove_modification,
                )
                mhcflurry_features = mhcflurry_generator.generate_features()
                psms.add_features(
                    mhcflurry_features,
                    psms_key=psms.peptide_column,
                    feature_key=mhcflurry_generator.id_column,
                    source="MHCflurry",
                )

                del mhcflurry_generator, mhcflurry_features
                gc.collect()

            elif generator_type == "NetMHCpan":
                #from optimhc.feature_generator.netMHCpan import NetMHCpanFeatureGenerator
                netmhcpan_generator = NetMHCpanFeatureGenerator(
                    unique_peptides,
                    alleles=allele,
                    mode=generator_params.get("mode", "best"),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    remove_modification=remove_modification,
                    n_processes=n_processes,
                    show_progress=show_progress,
                )
                netmhcpan_features = netmhcpan_generator.generate_features()
                psms.add_features(
                    netmhcpan_features,
                    psms_key=psms.peptide_column,
                    feature_key=netmhcpan_generator.id_column,
                    source="NetMHCpan",
                )

                del netmhcpan_generator, netmhcpan_features
                gc.collect()

            elif generator_type == "NetMHCIIpan":
                #from optimhc.feature_generator.netMHCIIpan import NetMHCIIpanFeatureGenerator
                netmhciipan_generator = NetMHCIIpanFeatureGenerator(
                    unique_peptides,
                    alleles=allele,
                    mode=generator_params.get("mode", "best"),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    remove_modification=remove_modification,
                    n_processes=n_processes,
                    show_progress=show_progress,
                )
                netmhciipan_features = netmhciipan_generator.generate_features()
                psms.add_features(
                    netmhciipan_features,
                    psms_key=psms.peptide_column,
                    feature_key=netmhciipan_generator.id_column,
                    source="NetMHCIIpan",
                )

                del netmhciipan_generator, netmhciipan_features
                gc.collect()

            elif generator_type == "DeepLC":
                from optimhc.feature_generator.DeepLC import DeepLCFeatureGenerator
                deeplc_generator = DeepLCFeatureGenerator(
                    psms,
                    calibration_criteria_column=generator_params.get(
                        "calibrationCriteria"
                    ),
                    lower_score_is_better=generator_params.get("lowerIsBetter"),
                    calibration_set_size=generator_params.get("calibrationSize", 0.1),
                    processes=n_processes, 
                    # TODO: Check here carefully
                    # Since DeepLC is GPU-based, 
                    # the processes here is not the same meaning as the n_processes in multi-threading
                    model_path=generator_params.get("model_path", None),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    mod_dict=mod_dict,
                )
                deeplc_features = deeplc_generator.generate_features()
                psms.add_features_by_index(
                    deeplc_features[deeplc_generator.feature_columns], source="DeepLC"
                )

                del deeplc_generator, deeplc_features
                gc.collect()

            elif generator_type == "SpectralSimilarity":
                from optimhc.feature_generator.spectral_similarity import SpectralSimilarityFeatureGenerator
                # Match PSMs with the spectra
                mzML_dir = generator_params.get("mzmlDir", None)
                if mzML_dir is None:
                    logger.error(
                        "mzML_dir is not provided for SpectralSimilarity feature generator."
                    )
                    continue

                pattern = generator_params.get("spectrumIdPattern", None)
                mz_file_names = []
                spectrum_ids = psms.spectrum_ids

                if pattern:
                    logger.info(
                        f"Using pattern: {pattern} to extract mzML file names from spectrum IDs."
                    )
                    for spectrum_id in spectrum_ids:
                        mz_file_names.append(re.match(pattern, spectrum_id).group(1))
                    logger.info(f"mzML file names: {list(set(mz_file_names))}")
                else:
                    logger.info("Spectrum ID pattern is not provided.")
                    if psms.ms_data_file_column is not None:
                        logger.info(
                            f"Trying to extract mzML file names from {psms.ms_data_file_column}"
                        )
                        logger.info(
                            f"MS data file format: {set(psms.psms[psms.ms_data_file_column])}"
                        )

                        for ms_data_file in psms.psms[psms.ms_data_file_column]:
                            mz_file_basename = os.path.basename(ms_data_file).split(
                                "."
                            )[0]
                            if mz_file_basename.endswith(".mzML"):
                                mz_file_basename = mz_file_basename[:-5]
                            elif mz_file_basename.endswith("mzML"):
                                mz_file_basename = mz_file_basename[:-4]
                            mz_file_names.append(mz_file_basename)

                        logger.info(f"mzML file names: {list(set(mz_file_names))}")
                    else:
                        logger.info("MS data file information is not provided.")
                        logger.info(
                            "Trying to use the default pattern: (.+?)\.\d+\.\d+\.\d+ to extract mzML file names from spectrum IDs."
                        )
                        for spectrum_id in spectrum_ids:
                            mz_file_names.append(
                                re.match("(.+?)\.\d+\.\d+\.\d+", spectrum_id).group(1)
                            )

                mz_file_paths = [
                    os.path.join(mzML_dir, f"{mz_file}.mzML")
                    for mz_file in mz_file_names
                ]
                mz_file_paths_set = set(mz_file_paths)
                logger.info(f"mz_file_paths: {mz_file_paths_set}")

                for mz_file_path in mz_file_paths_set:
                    if not os.path.exists(mz_file_path):
                        logger.error(f"mzML file not found: {mz_file_path}")
                        continue

                model_type = generator_params.get("model", None)
                if model_type is None:
                    logger.error(
                        "Model type is not provided for SpectralSimilarity feature generator."
                    )
                    raise ValueError(
                        "Model type is required for SpectralSimilarity feature generator."
                    )

                collision_energy = generator_params.get("collisionEnergy", None)
                instrument = generator_params.get("instrument", None)
                fragmentation_type = generator_params.get("fragmentationType", None)
                spectral_similarity_generator = SpectralSimilarityFeatureGenerator(
                    spectrum_ids=psms.spectrum_ids,
                    peptides=psms.peptides,
                    charges=psms.charges,
                    scan_ids=psms.scan_ids,
                    mz_file_paths=mz_file_paths,
                    model_type=generator_params.get("model"),
                    collision_energies=(
                        [collision_energy] * len(psms.peptides)
                        if collision_energy
                        else None
                    ),
                    instruments=(
                        [instrument] * len(psms.peptides) if instrument else None
                    ),
                    fragmentation_types=(
                        [fragmentation_type] * len(psms.peptides)
                        if fragmentation_type
                        else None
                    ),
                    remove_pre_nxt_aa=remove_pre_nxt_aa,
                    mod_dict=mod_dict,
                    url=generator_params.get("url"),
                    top_n=generator_params.get("numTopPeaks", 36),
                    tolerance_ppm=generator_params.get("tolerance", 20),
                )

                spectral_similarity_features = (
                    spectral_similarity_generator.generate_features()
                )
                psms.add_features(
                    spectral_similarity_features,
                    psms_key=[
                        psms.spectrum_column,
                        psms.peptide_column,
                        psms.charge_column,
                    ],
                    feature_key=spectral_similarity_generator.id_column,
                    source="SpectralSimilarity",
                )
                del (
                    spectral_similarity_generator,
                    spectral_similarity_features,
                    mz_file_paths,
                    mz_file_names,
                )
                gc.collect()

            else:
                logger.warning(f"Unknown feature generator: {generator_type}, skipping...")
