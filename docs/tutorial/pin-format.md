# PIN input support

PIN is a tab-separated interface rather than a fixed schema. Producers choose their own
feature names, charge encoding, spectrum identifiers, peptide notation, and protein layout.
OptiMHC detects columns from the header and reads one PIN file per call.

## Required fields

| Field | OptiMHC requirement |
| --- | --- |
| `SpecId` | Any string. It may be used to recover `run` and `rank`. |
| `Label` | Exactly `1` for target or `-1` for decoy. |
| `ScanNr` | Must be an integer. |
| `Peptide` | Must use one of the peptide forms below. |
| `Proteins` | Must be the last header field. Multiple proteins may use trailing tab fields or `;`. |
| Charge | Either one integer `Charge` column or supported one-hot columns. |

Optional metadata fields are matched case-insensitively:

- `rank`
- `filename` or `FileName`
- `ExpMass` and `CalcMass`
- `retention_time`, `retentiontime`, `ret_time`, or `rt`

`rank` is parsed separately but is always retained as a rescoring feature.

All other columns are considered possible rescoring features. Numeric columns are retained with
their original names. A column containing non-numeric values is omitted and reported as a
warning. NaN and infinite feature values are rejected later when Mokapot input is created.

## Retention time

Retention time is optional for reading and rescoring a PIN file, but DeepLC and FlashLFQ require
it. Comet does not normally write RT; Sage and MSFragger/FragPipe use `retentiontime`; Crux/Tide
may write `rt`. Generic Percolator input may use `rt` or `retentiontime`.

PIN does not define a universal RT unit. OptiMHC uses the largest value in the RT column to infer
the unit: values are treated as minutes and multiplied by 60 when the maximum is below 500;
otherwise they are treated as seconds. DeepLC and FlashLFQ receive the resulting seconds.

## Charge

The reader accepts:

- an integer `Charge` column;
- one-hot `Charge1`, `charge2`, or `charge_2` columns;
- one-hot columns such as `charge_7_or_more`, interpreted as charge 7.

Each row must select exactly one one-hot charge column. Sage's `z=2` through `z=6` columns and
numeric `z=other` field are not currently supported.

Input `Charge` or one-hot charge columns are also retained as rescoring features.

## Peptide and protein formats

Supported peptide examples include:

```text
K.PEPTIDE.R
K.PEPM[15.9949]K.R
K.n[42.0106]PEPTIDE.R
K.PEPTIDEc[-0.9840].R
PEPM[+15.9949]K
[+42.0106]-PEPTIDE
PEPTIDE-[-0.9840]
PEPM[UNIMOD:35]K
[UNIMOD:1]-PEPTIDE
PEPTIDE-[UNIMOD:2]
```

Flanking residues are optional. Numeric annotations are delta masses and must match a bundled
modification definition. The reader accepts Comet/MSFragger terminal markers (`n[mass]` and
`c[mass]`) and Sage terminal positions (`[mass]-` and `-[mass]`). Unimod annotations are retained
for OptiMHC PIN round trips. Crux symbol modifications are not supported.

Multiple proteins may be written either as trailing tab-separated values or as one
semicolon-separated `Proteins` value. OptiMHC stores both forms as a semicolon-separated string.

## Producer compatibility

| Producer | Current status | Relevant format |
| --- | --- | --- |
| MSFragger/FragPipe | Tested | Flanked delta-mass peptide, `ChargeN`, optional `retentiontime`. |
| OptiMHC PIN output | Tested round trip | Unimod peptide notation and integer `Charge`. |
| Comet | Format-compatible | Flanked delta-mass peptide, `Charge1...N`, trailing protein fields; RT is normally absent. |
| Crux/Tide | Partial | `ChargeN`, optional `rt`, and mass-bracket peptides work; symbol modifications do not. |
| OpenMS | Partial | Numeric peptide notation is supported; `ScanNr` must be an integer. |
| Sage | Partial | Numeric peptide notation is supported; string `ScanNr` and `z=...` charge fields are not. |

The optional Percolator `DefaultDirection` row is not supported. Compatibility is determined from
the actual header and values; OptiMHC does not require a fixed feature list for any producer.
