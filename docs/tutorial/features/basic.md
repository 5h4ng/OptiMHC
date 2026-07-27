# Basic Features

The **Basic** feature (`name: Basic`) computes simple sequence-level properties of each peptide. These features capture compositional and length characteristics that can help distinguish true identifications from random matches.

**Source name:** `Basic`

## Output Columns

| Column | Description |
|---|---|
| `length_diff_from_avg` | Peptide length minus the average peptide length across all PSMs |
| `abs_length_diff_from_avg` | Absolute value of the length difference |
| `unique_aa_count` | Number of distinct amino acid types in the peptide |
| `unique_aa_proportion` | Proportion of distinct amino acid types relative to peptide length |
| `shannon_entropy` | Shannon entropy of the amino acid frequency distribution |

## Computation

### Preprocessing

The feature uses the unmodified peptide sequence provided by the PSM reader.

### Length Features

Let \( L \) denote the length of a peptide and \( \bar{L} \) the average peptide length across the entire dataset of \( N \) PSMs:

\[
\bar{L} = \frac{1}{N} \sum_{i=1}^{N} L_i
\]

The length difference \( \delta \) and its absolute value are:

\[
\delta = L - \bar{L}, \quad |\delta| = |L - \bar{L}|
\]

### Unique Amino Acid Features

Let \( U \) denote the number of distinct amino acid types in a peptide. The unique amino acid proportion \( U_r \) is:

\[
U_r = \frac{U}{L}
\]

### Shannon Entropy

The Shannon entropy \( H \) quantifies the diversity of amino acid usage within a peptide:

\[
H = -\sum_{b \in \mathcal{A}} p_b \log_2(p_b)
\]

where \( \mathcal{A} \) is the set of unique amino acids in the peptide and \( p_b \) is the frequency of amino acid \( b \):

\[
p_b = \frac{n_b}{L}
\]

with \( n_b \) being the count of amino acid \( b \). A peptide composed of a single repeated amino acid has \( H = 0 \). A peptide with all unique amino acids has maximal entropy \( H = \log_2(L) \).

## Configuration

```yaml
featureGenerator:
  - name: Basic
```

No parameters are required.
