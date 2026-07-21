# Changelog

## Unreleased

### Breaking

- Replaced the configurable `PsmContainer` constructor with
  `PsmContainer(psm_df, feature_columns=...)`.
- Removed `spectrumIdPattern`. Spectral similarity now resolves each run as
  `<mzmlDir>/<run>.mzML`.

### Changed

- PIN and pepXML readers now produce the same internal PSM DataFrame schema.
- pepXML readers preserve every candidate rank.
- Multi-file inputs retain file order and use `run` plus `scan` for spectrum identity.
- Feature generators declare their output groups while existing `source` YAML remains supported.

### Fixed

- In-process Mokapot rescoring and exported PIN files now use the same projection.
- Charge one-hot features no longer cause an additional scalar charge feature.
- FlashLFQ output remains available through `toFlashLFQ`.
