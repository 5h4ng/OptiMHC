# Binding Affinity Score

OptiMHC integrates three external tools for predicting MHC binding affinity. These tools provide complementary estimates of peptide-MHC binding strength.

## Allele Names

Alleles are provided once in the top-level `allele` field and reused by binding affinity and PWM features:

```yaml
allele:
  - HLA-A*02:02
```

Use standard HLA names for NetMHCpan, NetMHCIIpan, and MHCflurry. Examples:

| MHC class | Example allele names |
|---|---|
| Class I | `HLA-A*02:01`, `HLA-A*02:02`, `HLA-B*07:02`, `HLA-C*07:02` |
| Class II DP/DQ | `HLA-DPA1*02:01-DPB1*01:01`, `HLA-DQA1*05:01-DQB1*02:01` |

NetMHCpan, NetMHCIIpan, and MHCflurry validate supported alleles through their own model packages or executables. PWM support is limited to alleles with matrices in `optimhc/PWMs/`; OptiMHC maps standard config names to PWM directory names internally. See [PWM Score](pwm.md#allele-support-and-format).

## NetMHCpan (MHC Class I)

**Config name:** `NetMHCpan` | **Source name:** `NetMHCpan`

[NetMHCpan 4.1](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) predicts peptide binding to MHC Class I molecules. OptiMHC uses NetMHCpan in binding affinity (BA) mode.

### Output Columns

| Column | Description |
|---|---|
| `netmhcpan_score` | Binding score for the selected allele |
| `netmhcpan_affinity` | Predicted binding affinity (nM) for the selected allele |
| `netmhcpan_percentile_rank` | Percentile rank for the selected allele |

### Computation

1. **Preprocess** peptide sequences: strip flanking amino acids and remove modifications.
2. **Filter** peptides to length 8–30 (the supported range for NetMHCpan).
3. **Predict** binding affinity for all peptide-allele combinations using `NetMHCpan 4.1 BA`.
4. **Allele selection**: for each peptide, select the allele \( a^* \) with the minimum percentile rank \( r \):

\[
a^* = \arg\min_{a} \; r(\text{peptide}, a)
\]

Missing values (peptides outside the length range) are filled with the column median.

### Configuration

```yaml
featureGenerator:
  - name: NetMHCpan
```

If `netMHCpan` is not on your `PATH`, provide the executable explicitly:

```yaml
featureGenerator:
  - name: NetMHCpan
    params:
      executablePath: /path/to/netMHCpan
```

!!! warning "External installation required"
    NetMHCpan is a standalone executable that must be downloaded from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/) and added to your `PATH`. See [Installation](../../getting-started/installation.md#netmhcpan-netmhciipan-setup) for details.

---

## NetMHCIIpan (MHC Class II)

**Config name:** `NetMHCIIpan` | **Source name:** `NetMHCIIpan`

[NetMHCIIpan 4.3](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) predicts peptide binding to MHC Class II molecules. It uses binding affinity (BA) mode.

### Output Columns

| Column | Description |
|---|---|
| `netmhciipan_score` | Binding score for the selected allele |
| `netmhciipan_affinity` | Predicted binding affinity (nM) for the selected allele |
| `netmhciipan_percentile_rank` | Percentile rank for the selected allele |

### Computation

1. **Preprocess** peptide sequences: strip flanking amino acids and remove modifications.
2. **Filter** peptides to length 9–50 (the supported range for NetMHCIIpan).
3. **Predict** binding for all peptide-allele combinations using `NetMHCIIpan 4.3 BA`.
4. **Allele selection**: for each peptide, select the allele with the minimum percentile rank.

Missing values are filled with the column median.

### Configuration

```yaml
featureGenerator:
  - name: NetMHCIIpan
```

If `netMHCIIpan` is not on your `PATH`, provide the executable explicitly:

```yaml
featureGenerator:
  - name: NetMHCIIpan
    params:
      executablePath: /path/to/netMHCIIpan
```

!!! warning "External installation required"
    NetMHCIIpan is a standalone executable that must be downloaded from [DTU Health Tech](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/) and added to your `PATH`. See [Installation](../../getting-started/installation.md#netmhcpan-netmhciipan-setup) for details.

---

## MHCflurry (MHC Class I)

**Config name:** `MHCflurry` | **Source name:** `MHCflurry`

[MHCflurry](https://github.com/openvax/mhcflurry) provides Class I MHC binding affinity predictions.

### Output Columns

| Column | Description |
|---|---|
| `mhcflurry_affinity` | Predicted binding affinity (nM) |
| `mhcflurry_processing_score` | Auxiliary MHCflurry model score |
| `mhcflurry_presentation_score` | Combined MHCflurry model score |
| `mhcflurry_presentation_percentile` | Percentile rank of the combined MHCflurry score |

### Computation

1. **Preprocess** peptide sequences: strip flanking amino acids and remove modifications.
2. **Filter** peptides to length 8–15 (the supported range for MHCflurry).
3. **Predict** MHCflurry outputs, including affinity and auxiliary model scores.

Missing values (peptides outside the 8–15 length range) are filled with the column median.

### Configuration

```yaml
featureGenerator:
  - name: MHCflurry
```

No additional parameters are required. The alleles are taken from the top-level `allele` configuration.

---

## Choosing Between Tools

| Tool | MHC Class | Peptide Length | Score | External Setup |
|---|---|---|---|---|
| NetMHCpan | I | 8–30 | Binding score, affinity, percentile rank | Manual download, add to `PATH` |
| NetMHCIIpan | II | 9–50 | Binding score, affinity, percentile rank | Manual download, add to `PATH` |
| MHCflurry | I | 8–15 | Affinity and auxiliary model scores | None (pip installed) |

For **MHC Class I** analyses, you can use both NetMHCpan and MHCflurry simultaneously — they provide complementary predictions. For **MHC Class II**, NetMHCIIpan is the only option among these three tools.
