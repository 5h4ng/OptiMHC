# Overlapping Peptide Score

The **OverlappingPeptide** feature (`name: OverlappingPeptide`) is designed for detecting **ladder-like peptide hotspots** — regions where related peptides form a characteristic nested/ladder pattern in the immunopeptidome.

The method is inspired by the **Overlap-Layout-Consensus (OLC)** algorithm from genome assembly. We adapt the first two stages — overlap detection and layout (contig construction) — to build an overlap graph of peptide sequences and assemble them into contigs. 

The core idea: if multiple peptides within the same contig are confidently identified, then other peptides in the same group are also likely to be genuine MHC ligands.

**Source name:** `OverlappingPeptide`

## Output Columns

| Column | Description |
|---|---|
| `contig_member_count` | Number of peptides in the contig that this peptide belongs to |
| `contig_extension_ratio` | How much the contig extends beyond the peptide itself |
| `contig_member_rank` | Rank of the contig by member count (largest contig = rank 1) |
| `contig_length` | Length of the assembled contig sequence |

## Computation

### Step 1: Preprocessing and Filtering

Base peptide sequences are filtered by length and entropy thresholds.

### Step 2: Redundancy Removal

Peptides that are exact substrings of longer peptides are identified. These **redundant** peptides are removed from the graph but are later mapped to the contig of their containing peptide. This prevents artificial inflation of contig membership.

### Step 3: Overlap Graph Construction

A directed graph \( G = (V, E) \) is built where:

- **Nodes** \( V \) represent non-redundant peptides.
- **Edge** \( A \to B \in E \) exists if a suffix of peptide \( A \) matches a prefix of peptide \( B \) with overlap length \( \geq \ell_\mathrm{min} \) (the `minOverlapLength` parameter).

For efficient overlap detection, a prefix index is built: for each peptide, all prefixes of length \( \ell_\mathrm{min} \) to \( |seq| \) are stored. Then for each peptide, its suffixes are looked up in the index to find overlapping partners.

When multiple overlaps exist between the same pair of peptides, only the **longest** overlap is kept. The edge weight stores the overlap length.

### Step 4: Transitive Edge Removal

To produce clean, non-branching paths, transitive edges are removed. If path \( A \to B \to C \) exists, the direct edge \( A \to C \) is removed (since the connection is already captured via \( B \)). This extends to longer transitive paths as well.

### Step 5: Contig Assembly

Contigs are identified as **non-branching paths** in the graph — chains of nodes where each internal node has exactly one predecessor and one successor.

Each contig is assembled into a consensus sequence by concatenating the non-overlapping portions. For a contig of peptides \( P_1, P_2, \ldots, P_k \) with overlap lengths \( o_1, o_2, \ldots, o_{k-1} \):

\[
C = P_1 \,\|\, P_2[o_1\!:] \,\|\, P_3[o_2\!:] \,\|\, \cdots \,\|\, P_k[o_{k-1}\!:]
\]

where \( P_i[o:] \) denotes the suffix of \( P_i \) beyond the first \( o \) characters, and \( \| \) denotes concatenation.

### Step 6: Peptide-to-Contig Mapping

Each peptide in a contig is assigned the contig's index. Redundant peptides (removed in Step 2) inherit the contig assignment of the peptide that contains them as a substring.

### Step 7: Feature Computation

Let \( n_c \) denote the number of peptides in a contig (including redundant ones), \( L_c \) the length of the assembled contig sequence, and \( L \) the length of the individual peptide.

The **contig member count** is simply \( n_c \).

The **contig length** is \( L_c \).

The **contig member rank** \( R_c \) is the rank of the contig by \( n_c \) in descending order (largest contig = rank 1, ties broken by minimum rank).

The **contig extension ratio** \( E_c \) measures how much the contig extends beyond the peptide:

\[
E_c = \frac{L_c - L}{L}
\]

A value of 0 means the peptide spans the entire contig; larger values indicate the contig extends significantly beyond the peptide.

### Missing Values

Peptides filtered out in Step 1 receive the median value of each generated
feature. This behavior is fixed when the generator is configured through YAML.

## Configuration

```yaml
featureGenerator:
  - name: OverlappingPeptide
    params:
      minOverlapLength: 7      # Minimum suffix-prefix overlap to form an edge
      minLength: 7             # Minimum peptide length to include
      maxLength: 20            # Maximum peptide length to include
      overlappingScore: expect # Search engine score column (for internal ranking)
```

| Parameter | Default | Description |
|---|---|---|
| `minOverlapLength` | `8` | Minimum overlap length for edge creation |
| `minLength` | `8` | Minimum peptide length to include in the graph |
| `maxLength` | `25` | Maximum peptide length to include |
| `overlappingScore` | Not set | Search score used to create contig-level aggregate features |

!!! tip
    For MHC Class I, typical settings are `minLength: 7`, `maxLength: 20`, `minOverlapLength: 7`. For MHC Class II, use wider ranges like `minLength: 9`, `maxLength: 50`, `minOverlapLength: 8` to accommodate longer peptides.
