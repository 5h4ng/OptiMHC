# Changelog

All notable user-visible changes to OptiMHC are documented here. The project currently
uses an `Unreleased` section while the package version remains `0.1.0`.

## [Unreleased]

### Breaking changes

- Replaced the configurable, source-aware `PsmContainer` with a canonical DataFrame-centered
  API. Construct it with `PsmContainer(df, feature_columns=...)`, access the table through
  `psms.df`, and attach features with `add_features(features, on=..., columns=...)`.
- Removed the mutable feature `source` registry from `PsmContainer`, together with the old
  constructor mappings, convenience projections, generic result mapping, index-based feature
  attachment, and container-owned PIN writer. Generator categories are instead recorded in a
  run-local feature manifest so existing experiment YAML remains compatible.
- Spectral similarity no longer uses `spectrumIdPattern`. Each canonical `run` resolves to
  `<mzmlDir>/<run>.mzML`, and `scan` selects the spectrum.

### CLI and configuration behavior

- The Click command surface is unchanged: `optimhc pipeline` and `optimhc experiment` remain
  the entry points, and `--inputFile` can still be repeated.
- Repeated `--inputFile` values and YAML `inputFile` lists are combined in caller-supplied
  order. Row order inside each file is preserved, then `psm_id` is reassigned sequentially.
- `rescore.seed` now defaults to integer `1` and is passed to both the selected model and
  `mokapot.brew`.
- Multi-file inputs must expose the same non-charge feature columns. Missing observed
  `charge_*` one-hot columns are the only schema difference filled automatically, using zero.
- Existing category-based experiment entries remain supported:

  ```yaml
  experiments:
    - name: baseline
      source: [Original, Basic]
      model: Percolator
  ```

- Each configured generator declares its output columns through the generator interface.
  `source` selects complete generator groups; `features` remains available as an advanced,
  mutually exclusive override for selecting individual DataFrame columns.
- `toFlashLFQ` remains supported and defaults to `true`. It writes accepted peptides to
  `<file_root>.FlashLFQ.txt` at `rescore.testFDR` through a standalone output adapter.

### Mokapot and exported PIN behavior

- In-process rescoring and exported PIN serialization now share one projection. Reading the
  exported PIN with Mokapot produces equivalent spectrum keys, rows, labels, peptides,
  proteins, and effective feature columns.
- Physical-spectrum competition uses `filename` plus `ScanNr`, corresponding to canonical
  `run` plus `scan`. Candidate `rank` and `SpecId` do not create additional spectra.
- Every candidate rank remains in the Mokapot dataset. `rank` is always a rescoring feature,
  even when an experiment requests a subset that omits it.
- When source PIN files already contain `charge_*` one-hot features, the shared projection omits
  the redundant scalar `Charge` column. Inputs without charge features still receive scalar
  `Charge`; internal rescoring and exported-PIN rescoring agree in both cases.
- Canonical retention time and experimental mass are intentionally omitted from the default
  Mokapot projection so they cannot change spectrum grouping. Calculated mass remains optional.
- Exported peptides use Unimod-accessioned HUPO-PSI ProForma annotations such as
  `M[UNIMOD:35]`. OptiMHC can read these exported PIN files back.
- With the same Mokapot version, model, FDR settings, and seed, the serialized PIN and
  in-process path now start from the same dataset. Mokapot 0.10's CLI accepts `--seed`, but does
  not pass it to `brew()`; exact reproducibility therefore requires the Python API with an
  explicit `rng` until that upstream CLI behavior changes.

### Reader and canonical-data behavior

- The required canonical columns are `psm_id`, `run`, `scan`, `rank`, `sequence`, `mods`,
  `mod_sites`, `charge`, `proteins`, and `is_decoy`.
- `read_pin()` and `read_pepxml()` read one file per call. Multi-file orchestration belongs to
  the pipeline.
- pepXML parsing preserves every `search_hit` in source order instead of retaining only the
  top-ranked candidate.
- PIN rank is taken from an explicit rank column, then from a trailing `_<rank>` in `SpecId`,
  and defaults to `1` only if neither is available.
- FragPipe/MSFragger PIN `retentiontime` values are converted from minutes to canonical seconds.
- PIN labels must be exactly `1` (target) or `-1` (decoy). Other values now fail instead of
  silently becoming targets.
- Declared rescoring features must be numeric and finite at the Mokapot boundary.
- Modifications are normalized with the vendored, attributed AlphaBase v1.9 modification table.
  Unknown masses or Unimod accessions fail explicitly.

### Fixed

- Fixed differing identification counts between OptiMHC rescoring and rescoring its exported
  PIN when the mismatch was caused by rank-dependent spectrum identities or different inferred
  feature columns.
- Fixed experiment-mode PIN files containing more features than the subset used in memory.
- Fixed duplicate peptide sequences breaking sequence-level feature generators when their
  predictions are identical.
- Fixed feature joins validating numeric conversions but storing the original object values.
- Fixed feature-importance plots receiving fewer column names than the fitted Mokapot model.
- Fixed exported MSBooster PIN files exposing both charge one-hot features and a redundant
  scalar `Charge`, which made Mokapot CLI train with one more feature than OptiMHC.
- Restored FlashLFQ output without returning export logic or Mokapot-private state to
  `PsmContainer`. Missing calculated masses are derived from canonical peptidoforms.

### Test coverage

- Contract tests cover the public `reader -> PsmContainer -> Mokapot projection` seam, including
  serialization round trips, ProForma round trips, invalid labels/features, feature joins, and
  multi-run/multi-rank competition.
- An actual seeded Mokapot regression verifies repeatable PSM and peptide confidence tables and no
  duplicate accepted `(filename, ScanNr)` values.
- Local integration tests cover MSFragger PIN and pepXML equivalence, MSBooster edited PIN files,
  Comet ranks 1-5, and two-file pipeline ordering.
- Large mzML fixtures remain gitignored. Run their integration check explicitly with:

  ```bash
  OPTIMHC_RUN_RAW_DATA=1 pytest -q tests/integration/test_real_psm_inputs.py -m raw_data
  ```

### Deferred

- Moving the legacy automatic p-value/E-value log transformation out of pepXML parsing remains
  deferred. Its numerical behavior is intentionally unchanged in this refactor.
