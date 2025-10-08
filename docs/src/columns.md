# Output Column Reference

Complete description of all columns in output files.

---

## Common Columns

These appear in multiple command outputs:

| Column | Type | Description |
|--------|------|-------------|
| `well` | Integer | Well identifier (1-based index in indices file) |
| `case` | String | Donor/case identifier |
| `name` | String | FASTQ read identifier |
| `genomic_sequence` | String | Full read sequence |
| `db_name` | String | Reference allele name (format: GENE*ALLELE, e.g., IGHV1-69*01) |
| `gene` | String | Gene name (extracted from db_name by splitting on "*") |
| `sequence` | String | Core allele sequence (may differ from db_seq for novel) |
| `count` | Integer | Number of supporting reads |
| `ratio` | Float | Allelic frequency (count / max_count in gene group) |
| `frequency` | Float | Alternate name for ratio in some commands |

---

## demultiplex Outputs

### {output}.tsv.gz

| Column | Description |
|--------|-------------|
| `well` | Well number from indices file (1-based) |
| `case` | Donor identifier from indices file |
| `name` | FASTQ read ID |
| `genomic_sequence` | Full read sequence |

### {output}.tsv.gz.log

| Column | Description |
|--------|-------------|
| `well`, `case` | Identifiers |
| `seq_count` | Number of reads in this well |
| `len_μ`, `len_σ` | Mean and standard deviation of read lengths |
| `len_q1`, `len_q2`, `len_q3` | Length quartiles (25%, 50%, 75%) |
| `len_min`, `len_max` | Length range |

---

## search exact Outputs

### Main Columns

| Column | Description |
|--------|-------------|
| `well`, `case`, `gene`, `db_name` | Identifiers |
| `sequence` | Core allele sequence |
| `count` | Reads matching this allele sequence (collapsed) |
| `full_count` | Reads matching sequence + flanks (uncollapsed) |
| `ratio` | Allelic ratio (count / max in gene) |
| `full_ratio` | Allelic ratio for full records |
| `flank_index` | Flank variant number (1 to --top) |

### Flank Columns (gene-dependent)

**V genes (3' RSS):**
- `prefix`: N bases before gene (--affix length, default 13)
- `heptamer`: 7bp RSS element (if --rss includes "heptamer")
- `spacer`: 23bp RSS spacer (if --rss includes "spacer")
- `nonamer`: 9bp RSS element (if --rss includes "nonamer")
- `suffix`: Extension sequence (if --extension specified)

**J genes (5' RSS):**
- `nonamer`, `spacer`, `heptamer`: 5' RSS elements
- `suffix`: N bases after gene (--affix length)
- `prefix`: Extension sequence (if --extension)

**D genes (both sides):**
- `pre_nonamer`, `pre_spacer`, `pre_heptamer`: 5' RSS (9+12+7 bp)
- `post_heptamer`, `post_spacer`, `post_nonamer`: 3' RSS (7+12+9 bp)
- Or `prefix`, `suffix` if --extension used

### Frequency Columns

| Column | Description |
|--------|-------------|
| `allele_freq` | Allele frequency within gene (per case): count / sum(counts in gene) |
| `allele_case_freq` | Allele frequency within case (all genes): count / case_count |
| `gene_case_freq` | Gene frequency within case: gene_count / case_count |
| `gene_count` | Total reads for this gene in this well/case |
| `case_count` | Total reads in this well/case |

### Cross-Case Comparison Columns

| Column | Description |
|--------|-------------|
| `cross_case_median_count` | Median count for this allele across all cases |
| `cross_case_median_allele_count` | Median allele count across cases |
| `cross_case_median_gene_count` | Median gene count across cases |
| `allele_to_cross_case_median_ratio` | count / cross_case_median_allele_count |
| `gene_to_cross_case_median_ratio` | gene_count / cross_case_median_gene_count |

**Use**: Identify donor-specific outliers (ratios <<1 suggest artifacts).

### Optional Columns

| Column | When Present | Description |
|--------|--------------|-------------|
| `isin_db` | --ref-fasta | "" (known) or "Novel" (not in reference) |
| `count_{refgene}_ratio` | --refgene | count / count of reference gene |
| `gene_count_{refgene}_ratio` | --refgene | gene_count / count of reference gene |

---

## search blast Outputs

### Core Columns

| Column | Description |
|--------|-------------|
| `qseqid` | Query sequence identifier |
| `sseqid` | Subject (database) sequence identifier |
| `gene` | Gene name (from sseqid) |
| `qseq` | Original query sequence |
| `db_seq` | Database reference sequence |
| `prefix`, `suffix` | Extension sequences |
| `aln_qseq` | Aligned query sequence (gaps removed, extensions trimmed) |
| `aln_mismatch` | Mismatches in aligned core (after trimming) |
| `corecov` | Core coverage = length(aln_qseq) / length(db_seq) |
| `isin_db` | Boolean: Is aln_qseq substring of known allele? |
| `allele_name` | Final name (original if exact match, or with "_S{hash} Novel") |
| `full_count` | Cluster size before collapsing |
| `full_frequency` | Frequency within well/case/gene |
| `count`, `frequency` | After aggregation |

### BLAST Standard Columns

| Column | Description |
|--------|-------------|
| `pident` | Percentage of identical matches |
| `nident` | Number of identical matches |
| `length` | Alignment length |
| `mismatch` | Number of mismatches |
| `gapopen` | Number of gap openings |
| `qcovs` | Query coverage per subject |
| `qcovhsp` | Query coverage per HSP |
| `qstart`, `qend` | Alignment start/end in query |
| `sstart`, `send` | Alignment start/end in subject |
| `qlen`, `slen` | Query and subject lengths |
| `evalue` | Expect value |
| `bitscore` | Bit score |
| `sstrand` | Subject strand |

---

## search hsmm Outputs

| Column | Description |
|--------|-------------|
| `well`, `case` | Identifiers |
| `sequence` | Detected D gene sequence |
| `count` | Reads with this D sequence (per well/case) |
| `ratio` | Allelic ratio within gene group (per well/case) |
| `gene` | Gene name (from nearest_db) |
| `allele_name` | Final allele name (reference or with _S{hash}) |

### RSS Elements

| Column | Length | Description |
|--------|--------|-------------|
| `pre_nonamer` | 9 bp | 5' nonamer |
| `pre_spacer` | 12 bp | 5' spacer |
| `pre_heptamer` | 7 bp | 5' heptamer |
| `post_heptamer` | 7 bp | 3' heptamer |
| `post_spacer` | 12 bp | 3' spacer |
| `post_nonamer` | 9 bp | 3' nonamer |

### Scoring Columns

| Column | Range | Description |
|--------|-------|-------------|
| `posterior_prob` | 0-1 | Confidence: exp(log_path_prob - log_total_prob) |
| `heptamer_prob_pre` | 0-1 | 5' heptamer PWM probability |
| `heptamer_prob_post` | 0-1 | 3' heptamer PWM probability |
| `heptamer_logp_pre` | ℝ | 5' heptamer log-probability |
| `heptamer_logp_post` | ℝ | 3' heptamer log-probability |
| `log_path_prob` | ℝ | Best path log-probability |
| `log_total_prob` | ℝ | Sum of all paths log-probability |

### Annotation Columns

| Column | Description |
|--------|-------------|
| `isin_db` | Boolean: Exact match to reference database? |
| `db_name` | Reference name if exact match, else "Novel" |
| `nearest_db` | Closest reference allele (by Levenshtein distance) |
| `nearest_db_dist` | Levenshtein distance to nearest_db |

---

## search heptamer Outputs

### Per-Read Output ({output}.tsv.gz)

| Column | Description |
|--------|-------------|
| `well`, `case`, `name`, `genomic_sequence` | From input |
| `db_name` | Reference allele name |
| `full_match` | Range of full database sequence in read (if found) |
| `trimmed_match` | Range of trimmed sequence in read |
| `db_length` | Length of database sequence |
| `full_length` | Length of extended sequence (trimmed + heptamer) |
| `full_sequence` | Extended sequence |
| `heptamer` | Identified heptamer sequence (7 bp) |

### Summary Output ({summary}.tsv)

| Column | Description |
|--------|-------------|
| `db_name` | Original reference name |
| `new_name` | Assigned name (with _S{hash} for novel) |
| `sequence` | Full extended sequence |
| `heptamer` | Heptamer sequence |
| `suffix` | Downstream sequence (includes heptamer) |
| `full_count` | Reads with exact match to untrimmed database sequence |
| `trimmed_count` | Reads with exact match to trimmed sequence |
| `allele_count` | Reads with this full_sequence + heptamer combination |
| `ratio` | allele_count / trimmed_count |
| `case` | Comma-separated list of cases |
| `full_length`, `db_length` | Sequence lengths |

---

## analyze association Outputs

### Edges Output ({edges}.tsv.gz)

| Column | Range | Description |
|--------|-------|-------------|
| `allele_a`, `allele_b` | - | Allele pair |
| `r` | -1 to +1 | Phi coefficient (Pearson correlation for binary data) |
| `r2` | 0 to 1 | Phi squared (≡ standardized LD r²) |
| `jaccard` | 0 to 1 | Jaccard index = n11/(n11+n10+n01) |
| `support` | Integer | n11 (donors with both alleles) |
| `similarity` | 0 to 1 | r or r2 after jaccard/support filtering |

### Clusters Output ({clusters}.tsv, optional)

| Column | Description |
|--------|-------------|
| `group_id` | Cluster number |
| `allele` | Allele name |
| `donors` | Comma-separated donor list |
| `n_donors` | Number of donors with this allele |
| `mean_n11`, `mean_n10`, `mean_n01`, `mean_n00` | Average contingency counts within cluster |
| `mean_r`, `max_r`, `min_r` | Phi coefficient statistics within cluster |

---

## analyze haplotype Outputs

| Column | Description |
|--------|-------------|
| `case` | Donor identifier |
| `gene` | Gene identifier |
| `genotype` | "homozygous", "heterozygous", or "uncertain" |
| `allele_1` | Primary allele (highest count) |
| `allele_2` | Secondary allele (2nd highest, empty if homozygous) |
| `count_1` | Count for allele_1 |
| `count_2` | Count for allele_2 (0 if homozygous) |
| `ratio` | Minor/major ratio = count_2 / count_1 |
| `total_count` | Total reads for this gene in this donor |
| `other_alleles` | Comma-separated additional alleles (if >2 present) |
| `novel_1`, `novel_2` | Boolean (if --novel-fasta provided) |

---

## analyze bwa Outputs

Added to input TSV (only rows mapping to target chromosome kept):

| Column | Description |
|--------|-------------|
| `position` | Chromosome and coordinate (format: "chr_name:pos") |
| `edit_distance` | NM field (number of mismatches in alignment) |
| `ref_sequence` | Reference genome sequence at position (orientation-adjusted) |
| `orientation` | "+" (forward strand) or "-" (reverse strand) |

All original input columns preserved.

---

## Column Interpretation Guide

### Count vs Frequency

- **count**: Absolute number of reads (integer)
- **ratio/frequency**: Relative proportion within group (0-1)
- **full_count**: Includes flanks (more specific, lower counts)
- **count**: Excludes flanks (collapsed, higher counts)

### Cross-Case Columns

Used to identify donor-specific artifacts:

```
allele_to_cross_case_median_ratio = count / median(count across all cases)
```

**Interpretation:**
- Ratio ≈ 1: Normal for this allele
- Ratio >> 1: Unusually high (possible artifact or true high expression)
- Ratio << 1: Unusually low (possible sequencing error specific to this donor)

### Phi Coefficient (r)

**Range:** -1 to +1

- **r = +1**: Perfect positive association (alleles always co-occur)
- **r = 0**: No association (independent)
- **r = -1**: Perfect negative association (mutually exclusive)

**Typical values:**
- r > 0.7: Strong positive association (likely same haplotype)
- r = 0.3-0.7: Moderate association
- r = 0-0.3: Weak association
- r < 0: Negative association (unlikely to co-occur)

### Posterior Probability (HSMM)

**Range:** 0 to 1

Confidence that the detected D segment is the best explanation:
```
posterior = exp(best_path_logprob - total_logprob)
```

**Interpretation:**
- >0.9: Very high confidence
- 0.7-0.9: High confidence (default threshold)
- 0.5-0.7: Moderate confidence
- <0.5: Low confidence (likely multiple competing explanations)

### Genotype Field

| Value | Meaning | Typical Ratio |
|-------|---------|---------------|
| homozygous | 1 allele OR minor allele below threshold | <0.1 |
| heterozygous | 2 alleles with ratio ≥ min-ratio | 0.1-0.7 |
| uncertain | >2 alleles present | varies |

**Note**: "uncertain" is unusual for diploid organisms and may indicate:
- Sequencing errors
- Gene duplication/triplication
- Mixed sample
- Low-quality data

---

## Data Type Reference

### String Columns
- Identifiers: `well`, `case`, `name`, `db_name`, `gene`, `allele_name`
- Sequences: `genomic_sequence`, `sequence`, `prefix`, `suffix`, RSS elements
- Annotations: `genotype`, `isin_db`, `orientation`

### Numeric Columns (Integer)
- Counts: `count`, `full_count`, `gene_count`, `case_count`, `support`
- Identifiers: `well`, `flank_index`, `group_id`, `n_donors`
- Distances: `aln_mismatch`, `nearest_db_dist`, `edit_distance`
- Lengths: `db_length`, `full_length`, `qlen`, `slen`

### Numeric Columns (Float)
- Frequencies: `ratio`, `full_ratio`, `frequency`, `allele_freq`, `*_case_freq`
- Correlations: `r`, `r2`, `jaccard`, `similarity`
- Probabilities: `posterior_prob`, `heptamer_prob_*`
- Coverage: `corecov`, `pident`, `qcovs`, `qcovhsp`
- Scores: `evalue`, `bitscore`
- Log-probs: `log_path_prob`, `log_total_prob`, `heptamer_logp_*`

### Boolean Columns
- `isin_db` (in some commands)
- `novel_1`, `novel_2` (if --novel-fasta)

**Note**: Some "boolean" columns use string values ("" vs "Novel") instead of true/false.

