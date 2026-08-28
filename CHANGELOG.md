# Changelog

## Unreleased

## 0.2.1 - 2026-08-29

### Fixed

- DeepLC calibration and retention-time prediction now run separately for each
  acquisition run, preventing pooled calibration across different chromatographic gradients.

## 0.2.0 - 2026-07-27

### Breaking

- Replaced the configurable `PsmContainer` constructor with
  `PsmContainer(psm_df, feature_columns=...)`.
- Removed `spectrumIdPattern`. Spectral similarity now resolves each run as
  `<mzmlDir>/<run>.mzML`.
- The legacy `modificationMap` key no longer controls PIN/pepXML parsing or DeepLC input.
  Readers now identify supported modifications using the bundled AlphaBase table before
  feature generation. Existing configurations are still accepted, but the key should be
  removed. DeepLC and spectral prediction now use reader-normalized modifications even when
  this key is absent.

### Changed

- PIN and pepXML readers now produce the same internal PSM DataFrame schema.
- PIN readers now ignore unknown non-numeric columns with a warning, recognize `rt` metadata,
  and convert retention times from minutes to seconds when the maximum value is below 500.
- PIN peptide parsing now accepts optional flanks and both Comet/MSFragger and Sage numeric
  terminal-modification syntax.
- pepXML readers preserve every candidate rank.
- Multi-file inputs retain file order and use `run` plus `scan` for spectrum identity.
- Feature generators declare their output groups while existing `source` YAML remains supported.

### Fixed

- In-process Mokapot rescoring and exported PIN files now use the same PIN DataFrame.
- Charge one-hot features no longer cause an additional scalar charge feature.
- FlashLFQ output remains available through `toFlashLFQ`.
