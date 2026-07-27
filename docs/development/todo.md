# Development TODO

This page lists planned improvements that are intentionally outside the current
PSM-container refactor.

## Read retention time from mzML

PIN files do not always provide retention time. When FlashLFQ output is enabled,
OptiMHC currently requires retention time in the PIN input.

Use mzML as a shared retention-time source so that downstream features and outputs
do not depend on search-engine-specific PIN columns. The input matching and
compatibility policy are still to be decided.

## Fix DeepLC retention-time calibration across runs

Tracked in
[issue #23](https://github.com/5h4ng/OptiMHC/issues/23).
The current DeepLC path calibrates all acquisition runs together:

- [DeepLC input preparation](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/feature/deeplc.py#L194-L208)
- [Global calibration and prediction](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/feature/deeplc.py#L244-L270)
- [Feature attachment by `psm_id`](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/feature/deeplc.py#L463-L467)

This can mix different chromatographic retention-time scales in a multi-run analysis.
Revisit the calibration design so that retention times remain comparable without
mixing run-specific behavior. The implementation and compatibility policy are
still to be decided.

## Simplify feature selection

The current feature-group compatibility layer is spread across the generator
interface, orchestration, pipeline state, and experiment configuration:

- [`BaseFeatureGenerator.feature_groups()`](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/feature/base_feature_generator.py#L39-L41)
- [Feature-group collection and validation](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/core/feature_generation.py#L105-L127)
- [Experiment feature selection](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/core/pipeline.py#L365-L370)
- [Legacy `source` configuration validation](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/core/config.py#L261-L275)

Revisit this API and reduce the amount of feature-selection state passed between
modules. The replacement and its migration path are still to be decided.

## Simplify PIN peptidoform handling

OptiMHC currently writes PIN peptides through
[`to_proforma()`](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/rescore/mokapot.py),
which produces Unimod annotations. The
[`read_pin()`](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/parser/pin.py)
path also contains dedicated support for these annotations.

Use numeric delta-mass annotations as the supported PIN modification format for
both reading and writing. OptiMHC-generated PIN files should use the Sage-style
peptidoform syntax. General ProForma and `UNIMOD:<id>` parsing do not need to be
part of the PIN reader. The compatibility impact and migration path are still
to be decided.

## Remove legacy peptide preprocessing options

PIN and pepXML readers now store the base peptide sequence separately from its
modifications. Remove the obsolete flank-removal configuration and duplicate
feature-level preprocessing after reviewing compatibility. If pepXML flanking
residues are needed in the future, store them as explicit PSM fields instead of
embedding them in the peptide sequence.

## Other follow-ups

- [ ] Separate search-score transformation from
  [`read_pepxml`](https://github.com/5h4ng/OptiMHC/blob/main/optimhc/parser/pepxml.py).
- [ ] Review the pinned Mokapot version when a suitable release is available.
