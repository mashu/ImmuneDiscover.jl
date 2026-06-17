# Parameter Guidelines

How to select appropriate thresholds and parameters for your analysis.

---

## Count Thresholds

Control minimum read support for alleles.

| Stringency | --mincount | Use Case |
|------------|------------|----------|
| **High confidence** | 10+ | Production genotyping, validated alleles |
| **Moderate** | 5 (default) | Balanced sensitivity/specificity |
| **Exploratory** | 1-3 | Rare alleles, pilot studies |

**Warning**: Count <5 increases false positive risk from sequencing errors.

---

## Frequency/Ratio Thresholds

Control allelic frequency within gene groups.

| Population Type | --min-allelic-ratio (exact) | --min-full-allelic-ratio / --min-peak-allelic-ratio (blast) | --min-ratio (haplotype) |
|-----------------|-----------------------------|-------------------------------------------------------------|-------------------------|
| **Homozygous-heavy** | 0.2 | 0.2 per-donor or peak | 0.2 |
| **Balanced** | 0.1 (default exact) | 0.1 per-donor; V preset uses peak 0.08 only | 0.1 (default) |
| **Low-expression** | 0.05 | 0.05 | 0.05 |

**Interpretation:**
- Ratio 0.5 = balanced heterozygous
- Ratio 0.1-0.3 = minor allele in heterozygous
- Ratio <0.1 = likely homozygous or artifact

---

## Distance Thresholds

Maximum mismatches for allele assignment.

| Gene | --max-blast-mismatch | --max-aln-mismatch | Rationale |
|------|---------------------|-------------------|-----------|
| **V genes** | 14 (`-g V`) | 14 (`-g V`) | BLAST cluster + trimmed-core distance for IGHV novel recovery |
| **D genes** | 20 (default) | 20 (default) | Highly variable, short |
| **J genes** | 10 (default) | 10 (default) | Moderately conserved |
| **Very strict** | 5 | 5 | High confidence only |

**Note**: `--max-blast-mismatch` filters the BLAST `mismatch` field (stored as `blast_mismatch` in the full table) before trimming; `--max-aln-mismatch` filters `core_aln_mismatch` in the output stage. They measure different things when affix extension is on.

---

## Extension Lengths

For BLAST and exact search.

| Gene | --forward | --reverse | Purpose |
|------|-----------|-----------|---------|
| **V genes** | 20 | 20 | RSS + stable affix alignment for extended DB (`-g V` preset) |
| **D genes** | 40 | 40 | Capture flanking V/J segments |
| **J genes** | 12 | 12 | Capture RSS heptamer |
| **No extension** | 0 | 0 | Use original sequences only |

**Warning**: Extensions <7 bp may cause alignment ambiguity.

---

## Coverage Thresholds

### Subject Coverage (--subjectcov)

Fraction of database sequence covered by alignment.

- **Default**: 0.1 (10%)
- **Strict**: 0.5+ (for full-length alleles)
- **D genes**: 0.25 (shorter, more variable)

### Core Coverage (--min-corecov)

Ratio of aligned sequence length to database sequence length after trimming.

- **Default**: 0.6 (60%)
- **V preset (`-g V`)**: 0.50 — keeps borderline full-length trims for novel alleles (raises false positives if lowered further)
- **Strict**: 0.8+ (require near-full coverage)
- **Permissive**: 0.4 (allow partial matches)

### Quality (--minquality)

Alignment quality = (length - mismatches) / length

- **Default**: 0.75 (75% identity after trimming)
- **V preset (`-g V`)**: 0.62 — slightly more permissive affix trimming for extended references
- **D genes**: 0.5 (more variable)
- **Strict**: 0.9+ (high identity only)

---

## Co-occurrence Analysis Parameters

### Allele inclusion

Control which alleles enter the co-occurrence network.

| Parameter | Conservative | Moderate | Exploratory |
|-----------|-------------|----------|-------------|
| `--min-donors` | 5+ | 2 (default) | 1 |

**Interpretation:**
- **min-donors**: Minimum number of distinct donors carrying the allele

All co-present allele pairs (support > 0) are written to the edges file with `rho`
(phi), `jaccard`, `support` (n11), `p_value` (hypergeometric enrichment) and the
Benjamini–Hochberg `q_value`; filter downstream on whichever column you need.

### Clustering

| Parameter | Value | Use Case |
|-----------|-------|----------|
| `--cluster-method complete` | default | Hierarchical clustering on rho (limits chaining) |
| `--cluster-method components` | | Connected components on rho ≥ threshold |
| `--cluster-method average/single` | | Hierarchical clustering on rho |
| `--cluster-threshold` | 0.7 (default) | Tight clusters (high confidence) |
| `--cluster-threshold` | 0.5 | Moderate clustering |
| `--cluster-threshold` | 0.3 | Loose clusters (exploratory) |
| `--min-cluster-size` | 3 (default) | Minimum cluster size to output |

---

## HSMM Parameters

### Training

- `--ratio` (default: 0.2): Balance between training data size and quality
  - Higher (0.3+): Fewer, higher-confidence training alleles
  - Lower (0.1): More training data, may include artifacts

- `--mincount` (default: 10): Ensures robust flank statistics
  - Minimum 5 recommended
  - Higher (20+) for cleaner PWM estimates

### Detection

- `--min-posterior` (default: 0.7): Detection confidence threshold
  - 0.9+: Very high confidence (may miss true positives)
  - 0.5-0.7: Moderate (balanced)
  - <0.5: Exploratory (more false positives)

- `--min-heptamer-prob-*` (default: 0.05): RSS quality filter
  - 0.1+: Strict RSS requirements
  - 0.01-0.05: Moderate
  - 0: Disable (rely only on posterior)

---

## Performance Tuning

### Multithreading

```bash
# Set before running commands
export JULIA_NUM_THREADS=16
```

**Commands that benefit:**
- `search exact`: Parallel read processing
- `discover hsmm`: Parallel detection scanning
- `discover blast`: BLAST uses threads internally
- `table` operations: Some use parallel processing

### Memory Management

**For large datasets:**
```bash
# Process plates individually
for demux in plate*_demux.tsv.gz; do
  immunediscover search exact $demux IGHV.fasta ${demux%.tsv.gz}_exact.tsv.gz -g V
done

# Then collect
immunediscover table collect "plate*_exact.tsv.gz" all_exact.tsv.gz
```

**For testing:**
```bash
# Limit reads
immunediscover search exact demux.tsv.gz IGHV.fasta test.tsv.gz -g V --limit 10000

# Reduce flank variants
immunediscover search exact demux.tsv.gz IGHV.fasta test.tsv.gz -g V --top 1
```

---

## Parameter Selection by Use Case

### Novel Allele Discovery (High Sensitivity)

```bash
immunediscover discover blast demux.tsv.gz IGHV.fasta blast.tsv.gz -g V \
  --min-fullcount 3 --min-peak-allelic-ratio 0.05 --minquality 0.7 --min-corecov 0.5
```

### Genotyping (High Specificity)

```bash
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv.gz -g V \
  --mincount 10 --minratio 0.1 --min-allele-cohort-fold 0.1
```

### Rare Allele Detection

```bash
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv.gz -g V \
  --mincount 1 --minratio 0.01 --min-allele-cohort-fold 0.01
```

### Haplotype Block Discovery

```bash
immunediscover analyze cooccurrence exact.tsv.gz \
  --min-donors 5 --cluster-method complete --cluster-threshold 0.7 --clusters clusters.tsv
```

---

## Common Parameter Combinations

### Conservative (Low False Positives)

```bash
# Exact search
-c 10 -f 0.1 --min-allele-cohort-fold 0.1 --min-gene-cohort-fold 0.1

# BLAST
-g V -c 10 -f 0.1 -d 5 -q 0.9 --min-corecov 0.8

# HSMM
--ratio 0.3 --mincount 20 --min-posterior 0.9 --out-mincount 20 --out-minratio 0.3
```

### Balanced (Default-like)

```bash
# Exact search
-c 5 -f 0.1 --min-allele-cohort-fold 0.05 --min-gene-cohort-fold 0.05

# BLAST
-g V  # Use preset defaults

# HSMM
--ratio 0.2 --mincount 10 --min-posterior 0.7 --out-mincount 10 --out-minratio 0.2
```

### Exploratory (High Sensitivity)

```bash
# Exact search
-c 1 -f 0.01 --min-allele-cohort-fold 0.01 --min-gene-cohort-fold 0.01

# BLAST
-g V -c 1 -f 0.01 -d 20 -q 0.5 --min-corecov 0.4

# HSMM
--ratio 0.1 --mincount 5 --min-posterior 0.5 --out-mincount 5 --out-minratio 0.1
```

