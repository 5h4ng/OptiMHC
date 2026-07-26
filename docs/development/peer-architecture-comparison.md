# OptiMHC PSM architecture: peer comparison and simplification decisions

Research date: 2026-07-17

## Scope and source snapshots

This note compares the current OptiMHC decisions and review findings with two primary-source snapshots:

- **Oktoberfest 0.11.1**, commit [`8ba32202`](https://github.com/wilhelm-lab/oktoberfest/tree/8ba32202a8b1231c124380a24a1070566b4114f0). The package version is declared in its [project metadata](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/pyproject.toml#L1-L16).
- **MS²Rescore 3.2.1**, commit [`d4121cad`](https://github.com/compomics/ms2rescore/tree/d4121cad220add428266e82596a23242fdee1784). Its metadata pins mokapot 0.10 and depends on `psm_utils` [in `pyproject.toml`](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/pyproject.toml#L1-L52).
- MS²Rescore delegates its internal PSM model and format adapters to **psm_utils 1.5.1**, commit [`c5960ee6`](https://github.com/compomics/psm_utils/tree/c5960ee61fe9fd6dd425965733ac669b9cb11ee8). This is an official CompOmics dependency, not a secondary source.
- Oktoberfest delegates PIN feature construction to **spectrum_fundamentals 0.10.0**, commit [`b9c8d1ec`](https://github.com/wilhelm-lab/spectrum_fundamentals/tree/b9c8d1ecdeb68cd9e41b5976e45c7f40b273192d), the version resolved by its lock file.

Only official repositories, their documentation in those repositories, and their tests were used. Statements that a feature is absent mean “not found in the inspected source tree at the pinned commit,” not a claim about unpublished branches.

## Executive result

The most useful peer pattern is **not** “rename the internal DataFrame to mokapot/PIN columns.” Neither peer does that:

- Oktoberfest owns an uppercase internal schema such as `RAW_FILE`, `SCAN_NUMBER`, `MODIFIED_SEQUENCE`, `PRECURSOR_CHARGE`, `SCORE`, and `REVERSE`; it converts that schema to `SpecId`, `Label`, `ScanNr`, `filename`, `Peptide`, and `Proteins` only when producing Percolator/mokapot input ([internal format](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/docs/internal_format.rst#L13-L42), [PIN mapping](https://github.com/wilhelm-lab/spectrum_fundamentals/blob/b9c8d1ecdeb68cd9e41b5976e45c7f40b273192d/spectrum_fundamentals/metrics/percolator.py#L341-L362)).
- MS²Rescore owns typed semantic fields (`spectrum_id`, `run`, `is_decoy`, `precursor_mz`, `retention_time`, `protein_list`, `rank`, `metadata`, `rescoring_features`) in `PSM`; its mokapot adapter creates and explicitly maps adapter-local columns ([PSM model](https://github.com/compomics/psm_utils/blob/c5960ee61fe9fd6dd425965733ac669b9cb11ee8/psm_utils/psm.py#L10-L34), [mokapot conversion](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/rescoring_engines/mokapot.py#L127-L186)).

**Recommendation:** keep a small canonical OptiMHC DataFrame schema, but make it an internal semantic schema. Put all `SpecId`/`Label`/`ScanNr` and mokapot optional-column naming in adapters. This adds one explicit seam but deletes pervasive column-role configuration. It is both simpler and more open to changing the rescoring engine.

## Observed peer architectures

### Oktoberfest

Oktoberfest's main in-memory module is `Spectra`, an `AnnData` subclass. PSM metadata live in `obs`; raw, predicted, and m/z fragment arrays live in shape-coupled sparse layers ([class and layer definitions](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/data/spectra.py#L42-L66), [column/layer mutation](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/data/spectra.py#L178-L258)). Search-engine readers convert into one internal DataFrame before spectrum matching. Feature construction is a comparatively large, fixed flow delegated to `spectrum_fundamentals.Percolator`; it writes PIN-like tab files, which are then read by mokapot or passed to Percolator ([generation flow](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/rescore/rescore.py#L22-L123), [mokapot call](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/rescore/rescore.py#L285-L331)).

This module has strong locality for fragment matrices, but the feature module is not very open to independently adding arbitrary generators: feature selection is expressed through flags and `additional_columns`, not a generator interface.

### MS²Rescore

MS²Rescore's main PSM module is `psm_utils.PSMList`, a collection of validated `PSM` records. Format readers are adapters into that model, and writers are adapters out. The orchestration module selects feature generators from a small explicit registry, asks each generator to mutate `rescoring_features`, checks completeness, converts to the selected rescoring engine, then writes typed results ([feature loop](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L51-L118), [generator interface and registry](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/feature_generators/base.py#L1-L30), [`FEATURE_GENERATORS`](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/feature_generators/__init__.py#L5-L20)).

The interface is more open to adding a generator than Oktoberfest's, although it is a source-code registry, not a third-party plugin system. The cost is a per-row object model and mutable nested dictionaries; OptiMHC need not copy that representation because its main structure is a DataFrame.

## Decision matrix

| Pending decision | Oktoberfest 0.11.1 | MS²Rescore 3.2.1 | Recommendation for OptiMHC | Simplicity / openness rationale |
|---|---|---|---|---|
| **Internal schema vs PIN/mokapot names** | Fixed internal uppercase schema; maps at PIN generation. | Typed semantic PSM model; adapter-local mokapot columns. | **Do not use PIN names as the internal schema.** Standardize a small semantic DataFrame schema and map once per adapter. | Deletes constructor column maps while keeping mokapot replaceable. |
| **Canonical required fields** | Broad workflow schema includes file, scan, sequence, charge, mass, score, decoy, length. | Only `peptidoform` and `spectrum_id` are structurally required; many fields are optional and workflow validation supplies constraints. | Require only fields the core path truly needs: stable `psm_id`, `spectrum_id`, `peptidoform`, `is_decoy`, and `protein_ids`; require charge/RT/run only at the generator/export seam that needs them. | Avoids false global requirements and optional-property crashes. |
| **PSM container representation** | `AnnData`: `obs` plus shape-coupled fragment layers. | `PSMList` of typed Pydantic records; recommends DataFrame conversion for vectorized work. | Keep DataFrame-first `PsmContainer`; do not adopt `AnnData` unless fragment matrices become core, and do not switch to per-row objects. | Preserves current performance and minimizes migration; a deep module can still own invariants. |
| **Parser adapters** | Search-engine-specific converters normalize before core processing. | `psm_utils.io` readers normalize many formats into `PSMList`; `parse_psms` then performs workflow validation ([read flow](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/parse_psms.py#L14-L47)). | PIN and PepXML parsers must return the same internal contract. Parser-specific column names must stop at the adapter. | Two adapters prove a real seam; shared contract tests provide leverage. |
| **PIN writer / round-trip** | PIN-like tab is an execution artifact, not the internal persistence format. | PIN writer emits required PIN columns plus selected features; it synthesizes `ScanNr` and does not preserve every internal field ([writer contract](https://github.com/compomics/psm_utils/blob/c5960ee61fe9fd6dd425965733ac669b9cb11ee8/psm_utils/io/percolator.py#L208-L280), [row mapping](https://github.com/compomics/psm_utils/blob/c5960ee61fe9fd6dd425965733ac669b9cb11ee8/psm_utils/io/percolator.py#L331-L357)). | Move `write_pin` to a PIN adapter. Test rescoring-semantic round-trip, not lossless internal-state round-trip. Use native TSV/Parquet for full persistence. | Removes format logic from the PSM module and avoids inflating PIN with unrelated metadata. |
| **Feature columns** | Generated columns live in a flat PIN table; `additional_columns` explicitly controls forwarded numeric input. | Each PSM has a flat `rescoring_features` dict; selected names become `feature:*` columns in the mokapot adapter. | Store one explicit `list[str]`/immutable set of feature columns; never infer features from “all non-metadata columns.” | Smallest interface and prevents RT/mass/rank becoming accidental features. |
| **Feature provenance / current `source`** | No general per-feature source state. Only fixed `original` and `rescore` artifacts distinguish bundles. | PSM `source` means file/search-engine provenance, not feature grouping. Generator→feature grouping is a run-local dict written to `.feature_names.tsv` and passed to reporting ([grouping](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L51-L60), [artifact](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L210-L216)). | **Delete feature `source` from `PsmContainer`.** If plots need grouping, return a run-level `FeatureManifest {generator: columns}` from orchestration; otherwise delete it entirely. | Removes mutable cross-cutting state without closing the generator extension seam. |
| **Ablation / experiment mode** | Has a fixed two-arm `original` vs prediction-enhanced comparison, run and plotted by the pipeline ([runner](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/runner.py#L590-L670)). No arbitrary generator-subset experiment module found. | No ablation framework found; configuration selects enabled generators. Report compares before/after rescoring. | **Retain existing category-based experiment YAML through a run-local feature manifest.** Keep exact column selection as an advanced alternative and keep all grouping state out of `PsmContainer`. | Preserves user workflows while keeping the canonical data module small. |
| **Feature-generator extension** | No general generator plugin seam; fixed `generate_features` flags and a large external `Percolator` feature implementation. | Small ABC (`feature_names`, `add_features`, `required_ms_data`) plus an explicit central registry. | Use one small protocol plus one explicit registry; no import side effects and no dynamic plugin framework yet. Prefer a generator to return a feature block rather than mutate container state. | Explicit registry is easy to read and change; returned blocks enable validation-before-commit. |
| **Metadata** | Metadata are ordinary `obs` columns; fragment arrays use dedicated layers. | Has free-form `metadata` and `provenance_data`, but core result fields are explicit. | Remove generic nested metadata from `PsmContainer`. Add a typed ordinary column only when two real consumers need it; keep generator-private intermediates inside the generator. | Better DataFrame copy semantics and locality; less speculative extensibility. |
| **Rescoring results** | Primarily file artifacts; results are not a generic mutable result registry on `Spectra`. | Explicit `score`, `qvalue`, `pep`, `rank`; peptide-level values are placed in metadata ([result update](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/rescoring_engines/mokapot.py#L218-L267)). | Delete generic `add_results`/result-column registry. Use explicit result columns returned by a rescorer adapter and validate one-to-one `psm_id`. | Explicit state is simpler and more testable; adapter remains replaceable. |
| **Atomic feature attachment** | No public generic join; AnnData layers must match observation shape. Feature tables are built positionally after resetting metadata index. | Generators mutate records, then core checks that every PSM has exactly the complete feature set and removes incomplete PSMs ([completeness check](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L93-L115)). No transaction guarantee was found. | Expose exactly one `with_features(block, on="psm_id")` operation: validate uniqueness, cardinality, finite numeric values, collisions, and complete match on a temporary DataFrame, then commit once. | Stronger than both peers and directly fixes partial mutation. |
| **Index alignment** | Resets metadata index before positional feature calculation; slices copy all AnnData axes together ([reset](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/rescore/rescore.py#L94-L103)). | Uses list position/artificial `index`, resets DataFrame index before mokapot, and maps confidence back by that artificial index ([conversion](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/rescoring_engines/mokapot.py#L142-L184), [mapping](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/rescoring_engines/mokapot.py#L218-L240)). | Never use ambient pandas index as identity. Create immutable unique `psm_id`; reset index after filtering only for presentation/performance. | Removes an entire class of silent alignment bugs while supporting reorder/filter operations. |
| **filename/run** | `RAW_FILE` is internal identity; mapped to PIN `filename`. | `run` is semantic optional state; mokapot adapter maps it via `filename_column`. | Choose internal name **`run`** and retire `ms_data_file`. Map it to mokapot `filename` only in the adapter. Require it only for multi-run competition or exports. | `run` matches the domain object, is shorter, and avoids confusing the file identity with an mzML path or open file handle. |
| **mass** | Internal `MASS`/`CALCULATED_MASS`; PIN `Mass` is a feature. Merged PIN `ExpMass` is deliberately repurposed as a per-(filename, ScanNr) group id, not experimental mass ([merge behavior](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/rescore/rescore.py#L197-L217)). | Computes `calcmass` from peptidoform and `expmass` from precursor m/z, then maps both explicitly. | Model `precursor_mz` (measured) and derive theoretical/experimental neutral mass in the mokapot adapter. Do not copy Oktoberfest's `ExpMass` workaround as domain semantics. | One measured source of truth; derived adapter data are easy to change. |
| **retention time** | `RETENTION_TIME` is internal metadata and becomes several generated numeric features (`RT`, `pred_RT`, `iRT`, `abs_rt_diff`) ([RT features](https://github.com/wilhelm-lab/spectrum_fundamentals/blob/b9c8d1ecdeb68cd9e41b5976e45c7f40b273192d/spectrum_fundamentals/metrics/percolator.py#L539-L559)). | `retention_time` is a semantic PSM field; DeepLC explicitly adds observed/predicted/difference features; mokapot receives RT as metadata. | RT is **metadata by default**, never automatically a rescoring feature. A generator may explicitly add named RT-derived feature columns. | Resolves PIN/PepXML disagreement and makes feature selection observable. |
| **charge** | Internal `PRECURSOR_CHARGE`; outputs one-hot `Charge1`…`Charge6` features. | Charge belongs to `peptidoform`; basic generator may add charge features; mokapot adapter explicitly passes `charge_column`. | Keep one semantic integer charge (or derive it from peptidoform). Charge one-hot columns are generated features, not core schema. | Avoids duplicated state while retaining modeling flexibility. |
| **Top-hit filtering** | Top-N selection exists for CE calibration, scoped to that analysis; general `Spectra` filters copy all axes ([CE selection](https://github.com/wilhelm-lab/oktoberfest/blob/8ba32202a8b1231c124380a24a1070566b4114f0/oktoberfest/predict/alignment.py#L12-L47)). | Rank is defined by collection/run/spectrum; input and output rank limits are workflow policy, and subsetting yields a new `PSMList` ([parse rank filter](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/parse_psms.py#L38-L42), [output filter](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L202-L207)). | Remove generic `get_top_hits` from the container unless the main path uses it. Put ranking in rescoring/analysis policy and preserve `psm_id`. | Keeps selection policy out of the data module; no gapped-index dependency. |
| **Validation and errors** | Property-based defaults plus targeted runtime checks; not one comprehensive schema. | JSON Schema/default cascade validates config; parsers raise contextual configuration errors and validate ID uniqueness/decoy presence ([config validation](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/config_parser.py#L120-L186), [PSM ID/decoy validation](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/parse_psms.py#L57-L79)). | Validate configuration once at entry; validate PSM invariants in `PsmContainer.__init__`; validate adapter-specific requirements at adapter entry. No “optional in constructor, indexed unconditionally” fields. | Errors stay local and defaults cannot drift from execution assumptions. |
| **Visualization/export** | Plotting is integrated after rescoring and many file artifacts are part of the runner. | HTML report and FlashLFQ are optional flags after core output; report errors are caught ([optional outputs](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/core.py#L178-L199)). | Make visualization and non-primary exports optional adapters invoked after the core result exists. They must consume result/manifest data, not mutate PSMs. | Core path becomes short; adapters stay independently removable/changeable. |
| **Test strategy** | Unit tests compare complete generated feature tables; integration tests exercise the runner. | Unit tests cover parsers, feature generators, rescoring adapters, reports, and `psm_utils` format round trips. | Add contract tests for every parser, atomic-attachment failure tests, filter/reorder identity tests, and one main-path integration test. Golden PIN test should assert only adapter contract. | The interface becomes the test surface; avoids testing private DataFrame mechanics. |

## The ten current correctness findings

The peers do not automatically solve all ten issues; the useful lesson is how to delete the interfaces that permit them.

| OptiMHC finding | Peer evidence / corresponding solution | Decision |
|---|---|---|
| 1. Top-hit filtering leaves gapped index; later index attachment silently creates `NaN`. | Oktoberfest resets before positional feature work; MS²Rescore resets and uses an artificial index, while ranking returns a new collection. | Immutable `psm_id`; eliminate index-based attachment; filter returns a new container or normalized copy. |
| 2. Failed duplicate-key feature attachment leaves rows, columns, and feature registration partially changed. | Neither peer exposes the same broad generic join. MS²Rescore checks completeness only after mutation, so it is not a transaction model to copy. | Validate a temporary merged DataFrame completely, then make one assignment; add a test that object state is byte-for-byte/equality unchanged after every failure. |
| 3. Returned DataFrame copies share nested dict/list objects; feature state is assigned by reference. | MS²Rescore intentionally uses nested record fields, but uses `deepcopy` before Percolator ID rewriting ([copy](https://github.com/compomics/ms2rescore/blob/d4121cad220add428266e82596a23242fdee1784/ms2rescore/rescoring_engines/percolator.py#L103-L116)). Oktoberfest keeps most tabular metadata scalar and matrices in layers. | Delete nested metadata/source structures from the main DataFrame. Return a deep copy only if object columns remain; prefer scalar/list-free core columns. |
| 4. Configured metadata column is ignored; merge cardinality/source detection unsafe. | Peers use named semantic fields or ordinary columns, not a configured “column containing metadata” with source detection. | Delete `metadata_column` and source detection rather than repair them. |
| 5. Normal same-named result join keys are rejected; result registry records/drops the key. | MS²Rescore has explicit result properties and maps back by stable artificial index/ID; Oktoberfest treats outputs as artifacts. | Delete generic result join/registry. Rescorer adapter returns exactly keyed result columns; container validates one-to-one and replaces/adds only allowed result names. |
| 6. Required schema columns can be deleted through feature mutation. | Peers distinguish semantic fields/metadata from feature fields. `psm_utils` PIN writer selects only requested `rescoring_features`. | Remove public arbitrary `drop_features`; if needed, allow only names in the explicit feature set and forbid core/result columns. |
| 7. Properties index optional columns even though constructor permits `None`. | `PSM` optional fields remain safe attributes; requirements are checked by the generator/adapter that needs them. | Constructor properties return `None`/absence safely; adapter uses `require_columns(...)` with a contextual error. Better: delete unused optional properties. |
| 8. Default configuration can raise `KeyError` because `featureGenerator` is absent. | MS²Rescore parses defaults and validation schema together before execution. | One typed/defaulted config model; tests must instantiate default config and run validation. Remove duplicated allowed-generator lists. |
| 9. Experiment `source` defaults to `None` but execution assumes iterable. | Neither peer has arbitrary source-subset experiments. | Resolve configured source names through the run-local generator manifest and fail explicitly on unknown categories. |
| 10. PIN and PepXML disagree whether RT is a feature. | Both peers separate semantic RT from format names; feature generation explicitly creates RT-derived features. | Canonical RT is metadata. Only explicitly registered RT-derived columns are features. Add parser contract assertion. |

## Edge-feature scope decisions

| Feature/scope | Oktoberfest | MS²Rescore | Recommendation |
|---|---|---|---|
| **MSBooster PIN style** | No MSBooster-specific writer found at the inspected commit. | No MSBooster-specific writer found; uses standard `psm_utils` Percolator writer. | Remove from core. Keep only if an active user/workflow exists, as a separate output adapter with its own golden test. |
| **FlashLFQ** | No FlashLFQ adapter found; quantification is a larger Picked Group FDR path. | Optional `write_flashlfq`; delegated to `psm_utils` writer, which maps run, sequence, theoretical mass, RT, charge, and proteins ([FlashLFQ mapping](https://github.com/compomics/psm_utils/blob/c5960ee61fe9fd6dd425965733ac669b9cb11ee8/psm_utils/io/flashlfq.py#L229-L245)). | Do not keep FlashLFQ logic in `PsmContainer`. Either delegate to a small adapter or delete until demanded. |
| **Binding-affinity (BA) Parquet intermediates** | No MHC binding-affinity feature found. Oktoberfest does use HDF5 as resumable workflow state. | No binding-affinity feature found. | Keep any BA cache entirely behind the BA generator seam. If it is only debugging residue, delete it; it must not shape the container interface. |
| **Visualization** | Rich plotting is integrated into the runner. | Optional report consumes PSMs plus generator→feature map. | Prefer MS²Rescore's data flow but make it a separate command/optional adapter if code reduction is the priority. Delete source-colored plots if source is deleted. |
| **Experiment mode** | Fixed original-vs-rescore comparison only. | No general experiment mode found. | Retain for backward compatibility, but resolve generator categories in orchestration; do not encode experiments or source groups in PSM state. |
| **Overlapping-peptide analysis** | No equivalent feature found. | No equivalent feature found. | If this is an OptiMHC differentiator, isolate it as one optional generator returning a flat feature block and no nested contig metadata. Otherwise remove it from the first simplified release. |

## Recommended target architecture

This is a recommendation, not an observation about the peers.

```text
PIN adapter ─┐
             ├─> PsmTable(DataFrame, feature_columns) ─> feature generators ─> rescorer adapter
PepXML adapter┘              │                                  │                    │
                             └─ invariants + atomic attach       └─ FeatureBlock      └─ ResultBlock

Optional, outside core: PIN export · FlashLFQ export · report · experiment runner
```

### Recommended internal contract

Use semantic names, not PIN names. The proposed first contract is:

- required identity/domain columns: `psm_id`, `spectrum_id`, `peptidoform`, `is_decoy`, `protein_ids`;
- conditionally required semantic columns: `run`, `scan_number`, `precursor_charge`, `precursor_mz`, `retention_time`, `search_score`;
- explicit state: `feature_columns: tuple[str, ...]`;
- explicit post-rescoring result columns: `score`, `q_value`, `pep`, `rank` (choose unambiguous pre/post names if both must coexist).

`psm_id` must be unique and immutable. It may be synthesized by an input adapter, but it must not be the pandas index contract. `scan` and `spectrum_id` should not both be globally required unless actual consumers prove both are necessary.

The mokapot/PIN adapter owns the external naming and value conversions:

| Internal semantic field | mokapot/PIN field or parameter | Rule |
|---|---|---|
| `psm_id` | `SpecId` or an adapter-local stable index | Preserve uniqueness; never use the pandas index as identity. |
| `run` | `filename` | Optional unless multi-run competition/export needs it. |
| `spectrum_id` | `spectrum_columns` / `SpecId` grouping input | The adapter decides whether multi-rank rows share a spectrum group. |
| `scan_number` | `ScanNr` / `scan_column` | Optional when a numeric scan exists; do not duplicate `spectrum_id` globally. |
| `is_decoy` | `Label` / `target_column` | Convert to PIN `-1/+1` or mokapot Boolean target semantics only here. |
| `peptidoform` | `Peptide` / `peptide_column` | Serialize modifications using the target format's convention. |
| `protein_ids` | `Proteins` / `protein_column` | Expand/join the internal protein collection only at export. |
| `precursor_charge` | `charge` / `charge_column` | Keep one integer internally; create one-hot charge features only when explicitly selected. |
| `retention_time` | `ret_time` / `rt_column` | Metadata, not a feature unless a generator creates explicit RT-derived columns. |
| derived masses | `calcmass`, `expmass` | Derive from `peptidoform`, `precursor_mz`, and charge in the adapter. |

### Recommended public interface

Keep the interface deliberately small:

```python
PsmTable.from_dataframe(df, feature_columns=...)
psms.dataframe()                    # defensive scalar-column copy/read view
psms.select(mask_or_ids)            # preserves psm_id
psms.with_features(feature_block)   # atomic, keyed by psm_id
psms.with_results(result_block)     # atomic, keyed by psm_id
psms.require_columns(names, context=...)
```

Do not retain generic `add_metadata`, `drop_source`, arbitrary `drop_features`, generic `add_results`, index-based feature attachment, or PIN serialization on the module.

## Final prioritized decisions

| Priority | Decision | Code simplification | Open to future change |
|---|---|---:|---:|
| P0 | Introduce immutable `psm_id`; ban pandas-index identity | Very high; eliminates two merge paths and many guards | High; all generators/adapters share one identity |
| P0 | Delete feature `source` from `PsmContainer` | Very high | High if an optional run-level manifest is retained |
| P0 | Delete experiment/ablation orchestration | Very high | Neutral: future runner can sit outside core |
| P0 | Make feature/result attachment atomic and keyed | Medium code addition, high defect deletion | High |
| P1 | Normalize PIN and PepXML into one semantic schema | High | High; new parsers become adapters |
| P1 | Move PIN output and mokapot mapping to adapters | High | Very high; rescorer/format can change |
| P1 | Replace feature dict-by-source with a flat explicit tuple/list | High | High |
| P1 | Delete nested generic metadata and generic result registry | High | Moderate; add typed columns only when demanded |
| P2 | Use a small explicit feature-generator protocol/registry | Medium | High without plugin-system complexity |
| P2 | Move/delete FlashLFQ, MSBooster, visualization, BA cache, overlapping-peptide extras according to actual users | Potentially very high | Use optional adapters/generators for the features retained |

The recommended bias is **deletion first, one deep PSM table module second, adapters third**. MS²Rescore demonstrates that format names and feature provenance do not need to live in the main PSM state; Oktoberfest demonstrates that a fixed internal schema can successfully feed both Percolator and mokapot. OptiMHC can take the simpler parts of each without adopting either peer's larger workflow surface.
