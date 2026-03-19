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

| Population Type | --minratio / --minfullratio | --min-ratio (haplotype) |
|-----------------|---------------------------|-------------------------|
| **Homozygous-heavy** | 0.2 | 0.2 |
| **Balanced** | 0.1 (default) | 0.1 (default) |
| **Low-expression** | 0.05 | 0.05 |

**Interpretation:**
- Ratio 0.5 = balanced heterozygous
- Ratio 0.1-0.3 = minor allele in heterozygous
- Ratio <0.1 = likely homozygous or artifact

---

## Distance Thresholds

Maximum mismatches for allele assignment.

| Gene | --maxdist | Rationale |
|------|-----------|-----------|
| **V genes** | 10 (default) | Moderately conserved |
| **D genes** | 20 (default) | Highly variable, short |
| **J genes** | 10 (default) | Moderately conserved |
| **Very strict** | 5 | High confidence only |

**Note**: Higher maxdist includes more divergent alleles but may capture pseudogenes.

---

## Extension Lengths

For BLAST and exact search.

| Gene | --forward | --reverse | Purpose |
|------|-----------|-----------|---------|
| **V genes** | 12 | 12 | Capture RSS heptamer |
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
- **Strict**: 0.8+ (require near-full coverage)
- **Permissive**: 0.4 (allow partial matches)

### Quality (--minquality)

Alignment quality = (length - mismatches) / length

- **Default**: 0.75 (75% identity after trimming)
- **D genes**: 0.5 (more variable)
- **Strict**: 0.9+ (high identity only)

---

## Association Analysis Parameters

### Support and Jaccard

Control edge inclusion in association network.

| Parameter | Conservative | Moderate | Exploratory |
|-----------|-------------|----------|-------------|
| `--min-support` (n11) | 10+ | 3 (default) | 1 |
| `--min-jaccard` | 0.5+ | 0.2 (default) | 0.1 |
| `--min-donors` | 5+ | 2 (default) | 1 |

**Interpretation:**
- **min-support**: Absolute co-occurrence count (n11)
- **min-jaccard**: Relative overlap = n11/(n11+n10+n01)
- **min-donors**: Total donors with allele (any context)

### Similarity and Clustering

| Parameter | Value | Use Case |
|-----------|-------|----------|
| `--similarity r` | default | Haplotype pairs (positive correlation only) |
| `--similarity r2` | | LD-style analysis (includes negative correlation) |
| `--threshold` | 0.5 (default) | Moderate clustering |
| `--threshold` | 0.7+ | Tight clusters (high confidence) |
| `--threshold` | 0.3 | Loose clusters (exploratory) |

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
  --minfullcount 3 --minfullratio 0.05 --minquality 0.7 --min-corecov 0.5
```

### Genotyping (High Specificity)

```bash
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv.gz -g V \
  --mincount 10 --minratio 0.1 --min-allele-mratio 0.1
```

### Rare Allele Detection

```bash
immunediscover search exact demux.tsv.gz IGHV.fasta exact.tsv.gz -g V \
  --mincount 1 --minratio 0.01 --min-allele-mratio 0.01
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
-c 10 -f 0.1 --min-allele-mratio 0.1 --min-gene-mratio 0.1

# BLAST
-g V -c 10 -f 0.1 -d 5 -q 0.9 --min-corecov 0.8

# HSMM
--ratio 0.3 --mincount 20 --min-posterior 0.9 --out-mincount 20 --out-minratio 0.3
```

### Balanced (Default-like)

```bash
# Exact search
-c 5 -f 0.1 --min-allele-mratio 0.05 --min-gene-mratio 0.05

# BLAST
-g V  # Use preset defaults

# HSMM
--ratio 0.2 --mincount 10 --min-posterior 0.7 --out-mincount 10 --out-minratio 0.2
```

### Exploratory (High Sensitivity)

```bash
# Exact search
-c 1 -f 0.01 --min-allele-mratio 0.01 --min-gene-mratio 0.01

# BLAST
-g V -c 1 -f 0.01 -d 20 -q 0.5 --min-corecov 0.4

# HSMM
--ratio 0.1 --mincount 5 --min-posterior 0.5 --out-mincount 5 --out-minratio 0.1
```

