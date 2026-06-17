# Command Reference

All immunediscover commands with complete parameter documentation.

---

## preprocess demultiplex

**Purpose:** Demultiplex indexed plate libraries into per-read metadata.

### Synopsis
```bash
immunediscover preprocess demultiplex <fastq> <indices> <output> [options]
```

### Arguments

**Required:**
- `fastq`: Input FASTQ file (gzipped or uncompressed)
- `indices`: TSV with `forward_index`, `reverse_index`, `case` columns (exactly 3, no others)
- `output`: Output TSV (auto-compressed to .gz)

**Optional:**
- `-l, --length` (default: 200): Minimum read length to keep
- `-s, --split`: Write per-case FASTQ files to `{fastq}_split/` directory
- `-f, --forwardarrayindex` (default: ""): Forward array index name for dual indexing
- `--case-filter-regex`: Regex to filter cases (e.g., `"[ACDERF]"`)

### Inputs/Outputs

**Input: indices.tsv**
```
forward_index	reverse_index	case
TGTTCGTGGT	AACATCCTCC	Donor1
GACACAAGGG	GACTCTCTAT	Donor2
```

**Output: {output}.tsv.gz**
- `well`: Well number (integer, 1-based)
- `case`: Donor identifier
- `name`: FASTQ read ID
- `genomic_sequence`: Full read sequence

**Output: {output}.tsv.gz.log**
- Statistics: `well`, `case`, `seq_count`, `len_μ`, `len_σ`, `len_q1/q2/q3`, `len_min/max`

### Examples

```bash
# Basic demultiplexing
immunediscover preprocess demultiplex plate1.fastq.gz indices.tsv demux.tsv.gz

# With length and case filtering
immunediscover preprocess demultiplex plate1.fastq.gz indices.tsv demux.tsv.gz \
  --length 250 --case-filter-regex "[ACDERF]"

# Batch processing
for f in *.fastq.gz; do
  base=${f%.fastq.gz}
  immunediscover preprocess demultiplex $f ${base}_indices.tsv ${base}_demux.tsv.gz
done
```

---

## search exact

**Purpose:** Exact match search against known alleles with RSS extraction and robust filtering.

### Synopsis
```bash
immunediscover search exact <tsv> <fasta> <output> -g <gene> [options]
```

### Arguments

**Required:**
- `tsv`: Demultiplexed TSV with `well`, `case`, `name`, `genomic_sequence`
- `fasta`: Reference allele database (FASTA format: GENE*ALLELE)
- `output`: Output TSV
- `-g, --gene`: Gene type - V, D, or J

**Count/Frequency Filters:**
- `-c, --mincount` (default: 5): Minimum read count
- `-f, --minratio` (default: 0.1): Minimum per-donor major ratio (`donor_major_ratio` and `donor_full_major_ratio`)
- `--min-allele-cohort-fold` (default: 0.05): Min fold of an allele's donor count vs its cohort median across donors
- `--min-gene-cohort-fold` (default: 0.05): Min fold of a gene's donor count vs its cohort median across donors
- `--min-recurrence` (default: 0, off): Require a candidate sequence in ≥ N donors (`n_donors`)
- `--min-seqlen` (default: 0, off): Drop candidates whose matched sequence is shorter than N nt
- `--min-peak-ratio` (default: 0.0, off): Require `peak_full_major_ratio` floor (best per-donor major ratio across donors). Complements `--minratio`.

All of these annotate rather than drop: every candidate is kept with `reject_reason` / `reject_stage`,
the run prints a per-filter kept/removed line at each stage, and two tables are written — the filtered
result and a full annotated table (`<output>.full.tsv.gz`).

**RSS Extraction:**
- `--rss` (default: "heptamer"): Extract RSS elements (heptamer, spacer, nonamer)
- `--extension`: Extension length in bp (replaces --rss)
- `-a, --affix` (default: 13): Bases from non-RSS side

**Advanced:**
- `-t, --top` (default: 1): Max flank variants per allele (1=collapsed mode)
- `-r, --refgene`: Reference gene(s) for ratio computation (space-separated)
- `--min-allelic-ratio` (default: 0.1): Primary filter on `donor_allelic_ratio` (allele fraction of gene reads per donor). Set 0 to disable.
- `-e, --expect`: Optional TSV overriding `--min-allelic-ratio` per gene/allele (columns: name, ratio); listed entries also skip `--minratio` dominance check
- `-d, --deletion`: TSV with gene_case_freq thresholds (columns: name, ratio)
- `--locus` (default: "IG"): Locus prefix for frequency calculations
- `--ref-fasta`: Reference FASTA to mark known vs novel (adds `isin_db` column)
- `--raw`: Path to save unfiltered results
- `-l, --limit` (default: 0): Limit input reads (0=unlimited, useful for testing)
- `-n, --noplot`: Disable unicode gene count boxplot

### Inputs/Outputs

**Input:** Demultiplexed TSV with `well`, `case`, `name`, `genomic_sequence`

**Output:** TSV with exact matches and flanks
- **Core**: `well`, `case`, `gene`, `db_name`, `sequence`
- **Counts**: `count`, `full_count`, `gene_count`, `case_count`, `n_reads_total`, `n_donors`, `flank_index`
- **Flanks** (schema follows the gene and mode searched — a V search carries no D columns):
  - RSS V: `prefix`, `sequence`, + selected `--rss` of `heptamer`/`spacer`/`nonamer`
  - RSS J: `suffix`, `sequence`, + selected `--rss` of `heptamer`/`spacer`/`nonamer`
  - RSS D: `pre_nonamer`/`pre_spacer`/`pre_heptamer` + `post_heptamer`/`post_spacer`/`post_nonamer`
  - Extension mode (any gene): `prefix`, `suffix`, `prefix_len`, `suffix_len`
- **Per-donor ratios**: `donor_major_ratio`, `donor_full_major_ratio`, `donor_allelic_ratio`, `peak_full_major_ratio`
- **Frequencies**: `gene_case_freq` (gene-usage / deletion), `*_cohort_fold`, `*_cohort_median`
- **Totals**: `gene_count`, `case_count`
- **Cohort fold-change**: `allele_cohort_median`, `gene_cohort_median`, `allele_cohort_fold`, `gene_cohort_fold` (count vs cohort-median across donors)
- **Other**: `chimera_score`, `reject_reason` / `reject_stage` (full table only)
- **Transparency**: `reject_reason`, `reject_stage` (full table only)
- **Markers**: `isin_db` (if --ref-fasta)
- **Ratios**: `count_{refgene}_ratio`, `gene_count_{refgene}_ratio` (if --refgene)

### Examples

```bash
# V gene search
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V

# D gene with extended flanks
immunediscover search exact demux.tsv.gz IGHD.fasta exact_D.tsv.gz -g D --extension 40

# Low stringency for rare alleles
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V \
  -c 1 -f 0.01 --min-allele-cohort-fold 0.01

# With reference ratios
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V \
  --refgene "IGHV1-69*01" "IGHV3-23*01"

# Mark novel sequences
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V \
  --ref-fasta IGHV_reference.fasta
```

### Notes

- **RSS extraction**: V genes have 3' RSS; J genes have 5' RSS; D genes have both sides
- **Filtering order**: Full records (sequence+flanks) filtered first, then collapsed records
- **Cross-case filtering**: Removes alleles abnormally low vs other donors (likely errors)
- **Expect/deletion files**: Override default thresholds for specific genes (e.g., pseudogenes)

---

## discover blast

**Purpose:** BLAST-based discovery with sequence extension, trimming, and identity clustering.

### Synopsis
```bash
immunediscover discover blast <input> <fasta> <output> -g <gene> [options]
```

### Arguments

**Required:**
- `input`: Demultiplexed TSV
- `fasta`: Database FASTA
- `output`: Output TSV
- `-g, --gene`: Gene preset (V, D, or J) - automatically sets optimal parameters

**Extension (primarily for D genes):**
- `--forward` (default: 20): Forward extension (bp)
- `--reverse` (default: 20): Reverse extension (bp)

**BLAST Settings:**
- `-a, --args` (default: gene-specific): Additional blastn arguments
- `-d, --max-blast-mismatch` (default: 20): Max BLAST `mismatch` per cluster (pre-trim hard filter)
- `--max-aln-mismatch` (default: 20): Max edit distance of trimmed core vs reference (output filter)
- `-e, --edge` (default: 0): Min nucleotides between gene and read end
- `-s, --subjectcov` (default: 0.1): Min subject coverage fraction

**Filtering:**
- `--min-fullcount` (default: 5): Minimum `full_count` per cluster
- `--min-count` (default: 0, off): Minimum collapsed `count` per (donor, allele, core)
- `--min-allelic-ratio` (default: 0, off): Min `allelic_ratio` (`count` ÷ max in donor+gene)
- `--min-full-allelic-ratio` (default: 0.1): Min `full_allelic_ratio` (`full_count` ÷ max in donor+gene). V preset sets 0 (off)
- `--min-peak-allelic-ratio` (default: 0, off): Min `peak_allelic_ratio` (max per-donor full ratio). V preset: 0.08
- `-l, --length` (default: 290): Min trimmed core length
- `-q, --minquality` (default: 0.75): Min alignment quality (1 - mismatch/length)
- `--min-corecov` (default: 0.6): Min ratio aligned_length/db_length
- `-i, --isin`: Keep truncated substrings of known alleles (default: filter out)
- `--keep-failed`: Keep rows with failed trimming (default: drop)

**Other:**
- `-p, --pseudo`: Pseudo-genes FASTA (prefixed with "P")
- `-o, --overwrite`: Overwrite BLAST cache files
- `-v, --verbose`: Save intermediate files
- `-G, --show-presets`: Display gene presets and exit

**Threading (not CLI flags):**
- **`blastn`** runs with `-num_threads` set from the environment variable **`BLAST_NUM_THREADS`** if set (positive integer), otherwise from the machine’s logical CPU count (`Sys.CPU_THREADS`).
- **Julia** uses its own thread pool for extension accumulation and trimming (`Folds`). Start the program with e.g. **`julia -t auto`** so that work is parallel; a single-threaded Julia process does not slow `blastn`, but post-BLAST steps will be slower.

### Gene Presets

**V Gene:**
- Extensions: forward=20, reverse=20
- Filtering: min-peak-allelic-ratio=0.08, min-full-allelic-ratio=0 (off), length=283, max-blast-mismatch=14, max-aln-mismatch=14, min-fullcount=5, minquality=0.62, min-corecov=0.50
- BLAST: `-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50`

**D Gene:**
- Extensions: forward=40, reverse=40
- Filtering: min-full-allelic-ratio=0.2, length=5, max-blast-mismatch=20, max-aln-mismatch=20, min-fullcount=10, edge=10, subjectcov=0.25, minquality=0.5
- BLAST: `-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5`

**J Gene:**
- Extensions: forward=12, reverse=12
- Filtering: length=10, max-blast-mismatch=10, max-aln-mismatch=10, min-fullcount=10
- BLAST: `-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10`

### Inputs/Outputs

**Input:** Demultiplexed TSV

**Outputs:**
- **{output}.tsv.gz**: Discovered alleles with columns for BLAST results, aligned sequences, coverage metrics, and novel allele names
- **`--work-dir`** (default: `.immunediscover` vs current working directory): all caches and intermediates — raw BLAST under `work-dir/blast/<16-hex>/hits.blast.gz` where `<16-hex>` hashes the **input TSV path + extended DB FASTA + `--args` + `--min-read-length`** (separate dirs per V/D input or DB/blastn settings; gene post-filters do not affect the key), plus combined/extended FASTA, `makeblastdb` files, and affix TSV
- **{fasta}-combined.fasta**: Combined database (pseudo + regular)
- **{fasta}-combined-extended.fasta**: Extended sequences (if extensions used)
- **{fasta}.affixes**: TSV with `name`, `prefix`, `suffix`

**Key columns**: `qseqid`, `sseqid`, `gene`, `qseq`, `db_seq`, `prefix`, `suffix`, `aln_qseq`, `blast_mismatch`, `core_aln_mismatch`, `corecov`, `isin_db`, `full_count`, `full_ratio`, `allele_name`

Plus standard BLAST columns: `pident`, `nident`, `length`, `mismatch`, `gapopen`, `qcovs`, `qcovhsp`, `qstart`, `qend`, `sstart`, `send`, `qlen`, `slen`, `evalue`, `bitscore`, `sstrand`

### Workflow

1. **Extension**: Extends database by most common flanks from reads (optional)
2. **BLAST**: Aligns reads to extended database; caches raw hits under `--work-dir` as streamed gzip (`work-dir/blast/<tag>/hits.blast.gz`)
3. **Re-alignment**: Trims extensions, re-aligns to original using Needleman-Wunsch
4. **Filtering**: Applies coverage, quality, frequency, distance filters
5. **Naming**: Assigns names with _S{hash} suffix for novel sequences

### Examples

```bash
# V gene with preset (caches under ./.immunediscover by default)
immunediscover discover blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz -g V

# D gene with preset
immunediscover discover blast demux_D.tsv.gz IGHD.fasta blast_D.tsv.gz -g D

# Override preset
immunediscover discover blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz \
  -g V --forward 20 --max-aln-mismatch 5

# With pseudo-genes
immunediscover discover blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz \
  -g V --pseudo IGHV_pseudo.fasta

# View presets
immunediscover discover blast --show-presets
```

### Notes

- **Extension purpose**: For D genes, captures flanking V/J to improve BLAST sensitivity
- **Cache**: BLAST results cached; use --overwrite to regenerate
- **Novel naming**: Format is `{closest_gene}_S{hash} Novel` where hash is MD5-based

---

## discover hsmm

**Purpose:** D gene detection using Hidden Semi-Markov Model trained on RSS flanks.

### Synopsis
```bash
immunediscover discover hsmm <tsv> <fasta> <output> [options]
```

### Arguments

**Required:**
- `tsv`: Demultiplexed TSV
- `fasta`: Known D allele database (for training)
- `output`: Output TSV

**Training Filters:**
- `-r, --ratio` (default: 0.2): Allelic ratio for training selection
- `-c, --mincount` (default: 10): Min count for training

**Model Parameters:**
- `--min-gene-len` (default: 10): Min D length for duration model (0=auto)
- `--max-gene-len` (default: 70): Max D length for duration model (0=auto)

**Detection Filters:**
- `--min-posterior` (default: 0.7): Min posterior probability of a sequence's best detection
- `--min-heptamer-prob-pre` (default: 0.05): Min 5' heptamer PWM probability (0=disable)
- `--min-heptamer-prob-post` (default: 0.05): Min 3' heptamer PWM probability (0=disable)

**Output Filters:**
- `--out-mincount` (default: 10): Min count for output
- `--out-minratio` (default: 0.2): Min allelic ratio for output

All of these annotate rather than drop. Detections are collapsed per `(well, case, sequence)` —
each cluster represented by its best (max-posterior) detection, with `count` = detections clearing
`--min-posterior` — then `--min-posterior` is applied as a "detection filter" stage and the count /
ratio / heptamer thresholds as an "output filter" stage. Every cluster (including low-posterior
ones) reaches the full table (`<output>.full.tsv.gz`) with `reject_reason` / `reject_stage`.

**Other:**
- `-l, --limit` (default: 0): Limit reads (0=unlimited)

### Inputs/Outputs

**Input:** Demultiplexed TSV with `well`, `case`, `name`, `genomic_sequence`

**Output: {output}.tsv.gz** - Detected D genes:
- **Identifiers**: `well`, `case`, `gene`, `allele_name`
- **Sequence**: `sequence`
- **RSS 5'**: `pre_nonamer` (9bp), `pre_spacer` (12bp), `pre_heptamer` (7bp)
- **RSS 3'**: `post_heptamer` (7bp), `post_spacer` (12bp), `post_nonamer` (9bp)
- **Scores**: `posterior_prob`, `heptamer_prob_pre`, `heptamer_prob_post`, `heptamer_logp_pre`, `heptamer_logp_post`, `log_path_prob`, `log_total_prob`
- **Annotation**: `isin_db`, `db_name`, `nearest_db`, `nearest_db_dist`
- **Metrics**: `count`, `ratio`

### Algorithm

The HSMM consists of 7 sequential states:
```
pre_nonamer(9) → pre_spacer(12) → pre_heptamer(7) → D_gene(variable) → 
post_heptamer(7) → post_spacer(12) → post_nonamer(9)
```

**Pipeline:**
1. Exact search known D alleles → extract complete RSS flanks
2. Per donor/gene: keep top 2 alleles with count≥mincount and ratio≥threshold
3. Train HSMM: PWMs for 6 RSS elements + IID emission + duration model for D gene
4. Scan all reads: find best D segment with flanks (Viterbi-like algorithm)
5. Filter by posterior probability, PWM probabilities, count, ratio

**Posterior probability** = exp(log_path_prob - log_total_prob) represents confidence in best detection vs all possible detections.

### Examples

```bash
# Basic D detection
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz

# Relaxed training
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --ratio 0.1 --mincount 5

# Strict detection
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --min-posterior 0.9 --min-heptamer-prob-pre 0.1 --min-heptamer-prob-post 0.1

# Custom gene length range
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --min-gene-len 8 --max-gene-len 50

# Test on subset
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz --limit 10000
```

### Notes

- Requires high-quality known D alleles across multiple donors for training
- Heptamer PWM filters assess RSS quality (0 disables filter)
- Detection phase is parallelized (use JULIA_NUM_THREADS)
- Nearest allele assigned by Levenshtein distance

---

## discover selftest

**Purpose:** Score how well discovery recovers known-novel alleles, and find which metric
threshold best separates true from false novel candidates.

### Synopsis
```bash
immunediscover discover selftest <discovery> <base> <truth> <output> [options]
```

### Arguments

**Required:**
- `discovery`: Discovery FULL table (`<output>.full.tsv.gz`) with `reject_reason` / `reject_stage`
- `base`: BASE reference FASTA used for discovery (known alleles only)
- `truth`: TRUTH FASTA (known + novel); truth-novel = sequences not in BASE
- `output`: TSV path for the per-allele recovery table

**Options:**
- `--seq-col` (default: `aln_qseq`): Discovery column holding the candidate core sequence
- `--no-substring`: Require exact sequence match (disable substring matching)
- `--metrics-output`: Optional TSV path to save the metric-separation table
- `-g, --gene` (optional): Gene preset used for `discover blast` — enables marginal filter shadowing audit at run thresholds

### Inputs/Outputs

**Input:** A discovery full table produced by running `discover blast` against the BASE reference,
plus the BASE and TRUTH FASTAs.

**Output: {output}** — per truth-novel allele: `allele`, `length`, `status`
(`recovered` / `rejected` / `missed`), `reject_stage`, and `reject_reason` (the exact filter label
from the full table, if rejected).

**Output (optional): {metrics-output}** — per metric: `metric`, `n_tp`, `n_fp`, `tp_median`,
`fp_median`, `direction` (`keep ≥` / `keep ≤`), `threshold`, `cli_suggestion` (matching
`discover blast` flag), `youden`, `tp_kept`, `fp_removed`.

### Algorithm

1. truth-novel = sequences in TRUTH not in BASE — the alleles discovery should find.
2. Classify each truth-novel allele as **recovered** (an accepted candidate matches), **rejected**
   (a candidate matched but was filtered — annotated with the stage), or **missed**.
3. Report **recall** (recovered / truth-novel) and **precision** (true / accepted-novel cores),
   plus outcome and rejected-by-stage bar plots.
4. Label every novel candidate row true/false and, per metric column, find the single threshold
   and keep-direction that best separates them by Youden's J (TP-rate − FP-rate). The most
   discriminative metric is reported first — that is the blast threshold to tune and to what value.

### Examples

```bash
# Run blast against the BASE reference to produce the full table, then score recovery
immunediscover discover blast demux.tsv.gz base.fasta disc.tsv.gz -g V
immunediscover discover selftest disc.tsv.gz.full.tsv.gz base.fasta truth.fasta recovery.tsv

# Save the metric-separation table for threshold tuning
immunediscover discover selftest disc.tsv.gz.full.tsv.gz base.fasta truth.fasta recovery.tsv \
  --metrics-output metrics.tsv
```

### Notes

- Consumes the FULL table, so it can show exactly which filter killed each true allele.
- Substring matching is on by default so a short truth allele matches a longer recovered core.

---

## search heptamer

**Purpose:** Identify heptamer RSS positions and extend V gene reads.

### Synopsis
```bash
immunediscover search heptamer <tsv> <fasta> <output> <summary> [options]
```

### Arguments

**Required:**
- `tsv`: Demultiplexed TSV
- `fasta`: V gene reference database
- `output`: Per-read output TSV
- `summary`: Collapsed summary TSV

**Optional:**
- `-j, --json` (default: "heptamers.json"): Heptamer dictionary JSON
- `-c, --chain` (default: "IGHV"): Chain (IGKV, IGLV, or IGHV)
- `-d, --maxdist` (default: 1): Max Hamming distance to heptamers
- `-b, --begin` (default: 0): Bases to trim from 5' start
- `-e, --end` (default: 8): Bases to trim from 3' end
- `-m, --mincount` (default: 1): Min count in summary
- `-r, --ratio` (default: 0.25): Min full/trimmed count ratio

### Heptamer Sequences

Auto-generated if heptamers.json missing:
- **IGHV**: CACAGTG, CACAATG, CACAGAG, CACGGTG, CACAGCG
- **IGKV**: CACAGTG, CACACTG, CACTGTG, CACGGTG, CACAATG, CACATTG
- **IGLV**: CACAGTG, CACGGTG, CATGGTG, CACGCTG, CACAGCG, CACAGTA, CATAGTG, CACAATG

### Inputs/Outputs

**Output: {output}.tsv.gz** - Per-read:
- `well`, `case`, `name`, `genomic_sequence`, `db_name`
- `full_match`, `trimmed_match`: Match ranges
- `db_length`, `full_length`: Lengths
- `full_sequence`, `heptamer`: Extended sequence and identified heptamer

**Output: {summary}.tsv** - Collapsed:
- `db_name`, `new_name`, `sequence`, `heptamer`, `suffix`
- `full_count`, `trimmed_count`, `allele_count`, `ratio`
- `case`: Comma-separated list
- `full_length`, `db_length`

### Examples

```bash
# Basic IGHV search
immunediscover search heptamer demux.tsv.gz IGHV.fasta hept.tsv.gz summary.tsv

# Allow 2 mismatches
immunediscover search heptamer demux.tsv.gz IGHV.fasta hept.tsv.gz summary.tsv \
  --maxdist 2

# Strict summary filtering
immunediscover search heptamer demux.tsv.gz IGHV.fasta hept.tsv.gz summary.tsv \
  --mincount 10 --ratio 0.5
```

---

## analyze cooccurrence

**Purpose:** Compute pairwise allele co-occurrence and associations across donors.

### Methodology Note

This analysis computes **phi coefficient** (correlation-based association) rather than **classical linkage disequilibrium (LD)** because:

**Our data is unphased**: We know which alleles are present in each donor, but not which chromosome (maternal vs. paternal) they're on. Classical LD measures like D (= pAB - pA×pB) and D' require phased haplotype data.

**What we compute:**
- **Phi coefficient (φ)**: Pearson's correlation for binary presence/absence across donors
- **φ²** (r²): Mathematically equivalent to standardized LD r² used in population genetics
- Measures: "Do alleles A and B co-occur more than expected by chance?"

**Note**: φ² is mathematically equivalent to the standardized LD r² used in population genetics, so the reported `rho` (φ) gives LD-comparable information from unphased genotype data.

### Synopsis
```bash
immunediscover analyze cooccurrence <input> [options]
```

### Arguments

**Required:**
- `input`: TSV with case and allele columns (the edges file is written automatically next to it)

**Column Names:**
- `-C, --case-col` (default: "case"): Donor column
- `-A, --allele-col` (default: "db_name"): Allele column

**Filtering:**
- `-m, --min-donors` (default: 2): Minimum distinct donors required to include an allele

**Clustering:**
- `--cluster-method` (default: "complete"): `complete`, `components`, `average`, or `single`
- `--cluster-threshold` (default: 0.7): Similarity (rho) cutoff for clustering
- `--min-cluster-size` (default: 3): Minimum cluster size to output
- `--clusters`: Optional path to write clusters TSV
- `--debug-triangles`: Print rho-based triangle diagnostics at the cluster threshold

### Inputs/Outputs

**Input:** Any TSV with at least the case and allele columns (defaults: `case`, `db_name`)

**Output: {input}_edges.tsv** - one row per co-present allele pair (support > 0):
- `allele_a`, `allele_b`: Allele pair
- `rho`: Phi coefficient (-1 to +1)
- `jaccard`: Jaccard index (0 to 1)
- `support`: n11 (donors carrying both)
- `p_value`: Hypergeometric enrichment p-value
- `q_value`: Benjamini–Hochberg adjusted p-value

**Output: {input}_clusters.tsv** — partial-haplotype blocks (always written, sorted by `group_id`):
- `group_id`: Block id (`0` = unclustered)
- `group_size`: Alleles in the block
- `allele`, `donors`, `n_donors`

Override path with `--clusters`.

### Metrics

**Contingency Table (alleles A, B):**
- n11: Donors with both
- n10: Donors with A only
- n01: Donors with B only
- n00: Donors with neither

**Phi Coefficient:**
```
φ = (n11×n00 - n10×n01) / sqrt((n11+n10)(n11+n01)(n10+n00)(n01+n00))
```
Pearson's r for binary data. Range: -1 to +1.

**Jaccard Index:**
```
J = n11 / (n11 + n10 + n01)
```
Overlap ignoring double-absence. Range: 0 to 1.

φ² is equivalent to standardized LD r² but computed from unphased data.

### Examples

```bash
# Basic
immunediscover analyze cooccurrence exact_V.tsv.gz

# Strict filtering
immunediscover analyze cooccurrence exact_V.tsv.gz \
  --min-donors 5

# With clustering
immunediscover analyze cooccurrence exact_V.tsv.gz \
  --cluster-method complete --cluster-threshold 0.7 --clusters clusters.tsv

# Custom columns
immunediscover analyze cooccurrence input.tsv.gz \
  --case-col donor_id --allele-col allele_name
```

### Notes

- **r vs r²**: "r" keeps only positive correlations; "r²" includes negative
- **Unphased data**: Measures population co-occurrence without chromosome phase info
- **Clustering**: Complete-linkage hierarchical on distance = 1 - similarity
- **Use cases**: Haplotype blocks, novel allele validation, artifact detection

---

## analyze haplotype

**Purpose:** Infer homozygous vs heterozygous genotypes using diploid assumptions.

### Synopsis
```bash
immunediscover analyze haplotype <input> <output> [options]
```

### Arguments

**Required:**
- `input`: TSV with case, gene, allele, count columns
- `output`: Output TSV

**Column Names:**
- `-C, --case-col` (default: "case")
- `-A, --allele-col` (default: "db_name")
- `-G, --gene-col` (default: "gene")

**Filtering:**
- `-c, --mincount` (default: 5): Min count threshold
- `-r, --min-ratio` (default: 0.1): Min minor/major allele ratio for heterozygous

**Optional:**
- `-f, --novel-fasta`: Novel alleles FASTA (adds `novel_1`, `novel_2` columns)

### Inputs/Outputs

**Input:** TSV with `case`, `gene`, `db_name` (allele), `count` columns. Typically from `search exact` or `discover blast`.

**Output: {output}.tsv**:
- `case`, `gene`, `genotype` (homozygous/heterozygous/duplication)
- `allele_1`, `allele_2`: Primary and secondary alleles
- `count_1`, `count_2`: Respective counts
- `ratio`: Minor/major ratio (count_2 / count_1)
- `total_count`: Total for this gene in donor
- `other_alleles`: Additional alleles if >2 present
- `novel_1`, `novel_2`: Boolean (if --novel-fasta)

**Genotype Classification:**
- **Homozygous**: 1 allele OR ratio < min-ratio
- **Heterozygous**: 2 alleles with ratio ≥ min-ratio
- **Uncertain**: >2 alleles (unusual for diploid)

### Examples

```bash
# Basic
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv

# Mark novel alleles
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv \
  --novel-fasta novel_V.fasta

# Relaxed threshold
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv \
  --min-ratio 0.05 --mincount 3

# Custom columns
immunediscover analyze haplotype input.tsv.gz haplotypes.tsv \
  --case-col donor --allele-col allele_name --gene-col gene_id
```

### Notes

- Assumes diploid genetics (≤2 functional alleles per gene)
- Ratio <0.1 typically indicates homozygous or sequencing artifact
- Ratio ~0.5 indicates balanced heterozygous genotype

---

## search bwa

**Purpose:** BWA genome mapping QC to retain chromosome-specific sequences.

### Synopsis
```bash
immunediscover search bwa <tsv> <output> <genome...> [options]
```

### Arguments

**Required:**
- `tsv`: Input TSV with allele columns
- `output`: Filtered output TSV
- `genome`: BWA-indexed genome FASTA(s) - can provide multiple files

**Optional:**
- `-c, --chromosome` (default: "chromosome 14"): Target chromosome string
- `-t, --tag` (default: "(.*Primary Assembly.*)|(.*alternate locus.*)"): Chromosome description regex
- `-n, --colname` (default: "best_name"): Allele name column
- `-s, --colseq` (default: ["prefix", "best_aln", "suffix"]): Sequence columns to concatenate

### Inputs/Outputs

**Input:** TSV with allele name and sequence column(s)

**Setup:** Index genome once:
```bash
bwa index genome.fasta
```

**Output: {output}.tsv.gz** - Input rows mapped to target chromosome with added:
- `position`: "chromosome:coordinate"
- `edit_distance`: NM (mismatches)
- `ref_sequence`: Reference genome sequence (orientation-adjusted)
- `orientation`: "+" (forward) or "-" (reverse)

**Side output: /tmp/discarded.tsv** - Off-target allele names

### Examples

```bash
# Human chromosome 14
bwa index GCF_000001405.25.fasta
immunediscover search bwa candidates.tsv.gz filtered.tsv.gz \
  GCF_000001405.25.fasta --chromosome "chromosome 14"

# Mouse chromosome 12
immunediscover search bwa candidates.tsv.gz filtered.tsv.gz \
  mouse_genome.fasta --chromosome "chromosome 12"

# Custom columns
immunediscover search bwa candidates.tsv.gz filtered.tsv.gz \
  genome.fasta --colname allele_id --colseq sequence
```

### Notes

- Genome must be BWA-indexed before use
- Can provide multiple FASTA files (e.g., separate chromosomes)
- `ref_sequence` is reverse-complemented if alignment on reverse strand
- Essential QC step before validating novel alleles

---

## table Commands

### table outerjoin / leftjoin

**Purpose:** Join two TSV files on key columns.

**Synopsis:**
```bash
immunediscover table {outerjoin|leftjoin} <left> <right> <output> --keys <cols> [options]
```

**Arguments:**
- `left`, `right`, `output`: File paths
- `-k, --keys`: Comma-separated join columns
- `--left-keys`, `--right-keys`: Different keys per file (optional)
- `--left-prefix`, `--right-prefix`: Column prefixes (optional)
- `--left-select`, `--right-select`: Column subsets (optional)

**Example:**
```bash
immunediscover table outerjoin exact.tsv.gz haplotypes.tsv merged.tsv.gz \
  --keys case,gene --left-prefix exact --right-prefix haplo
```

---

### table transform

**Purpose:** Regex-based column transformation with capture groups.

**Synopsis:**
```bash
immunediscover table transform <input> <output> -c <col> -p <pattern> -r <replacement>
```

**Arguments:**
- `-c, --column`: Target column
- `-p, --pattern`: Regex with capture groups
- `-r, --replacement`: Replacement using \1, \2, etc.
- `--new-column`: New column name (preserves original)

**Example:**
```bash
immunediscover table transform input.tsv.gz output.tsv.gz \
  --column db_name --pattern 'IGHV(.*)\\*(.*)' --replacement 'V\\1-\\2'
```

---

### table aggregate

**Purpose:** Group by columns and count.

**Synopsis:**
```bash
immunediscover table aggregate <input> <output> --group-by <cols> [options]
```

**Arguments:**
- `-g, --group-by`: Comma-separated grouping columns
- `-k, --keep-columns`: Additional columns to keep (first value)
- `-c, --count-column` (default: "count"): Count column name

---

### table unique

**Purpose:** Remove duplicate rows by specified columns.

**Synopsis:**
```bash
immunediscover table unique <input> <output> --columns <cols>
```

---

### table sort

**Purpose:** Sort by columns.

**Synopsis:**
```bash
immunediscover table sort <input> <output> --columns <cols> [--reverse]
```

**Arguments:**
- `-c, --columns`: Comma-separated sort columns (priority order)
- `-r, --reverse`: Descending order (default: ascending)

---

### table filter

**Purpose:** Filter rows by regex pattern or numeric threshold.

**Synopsis:**
```bash
# Regex mode
immunediscover table filter <input> <output> -c <col> --pattern <regex>

# Numeric mode
immunediscover table filter <input> <output> -c <col> --operator <op> --threshold <val>
```

**Arguments:**
- `-c, --column`: Target column
- `--pattern`: Regex pattern (string filtering)
- `--operator`: Comparison (<, <=, >=, >)
- `--threshold`: Numeric value

**Examples:**
```bash
# Keep Novel alleles
immunediscover table filter blast.tsv.gz novel.tsv.gz \
  --column allele_name --pattern "Novel"

# Keep high-count
immunediscover table filter blast.tsv.gz filtered.tsv.gz \
  --column count --operator ">=" --threshold 10
```

---

### table select

**Purpose:** Project subset of columns.

**Synopsis:**
```bash
immunediscover table select <input> <output> --columns <cols>
```

**Example:**
```bash
immunediscover table select blast.tsv.gz simple.tsv.gz \
  --columns case,gene,allele_name,count
```

---

### table fasta

**Purpose:** Export TSV sequences to FASTA format.

**Synopsis:**
```bash
immunediscover table fasta <input> <output> [options]
```

**Arguments:**
- `-n, --colname` (default: "allele_name"): Name column
- `-s, --colseq` (default: "seq"): Sequence column
- `-d, --coldesc`: Description column (optional)
- `-f, --filter`: Regex to filter names
- `-c, --cleanup`: Regex to remove from names (e.g., " Novel")
- `--desc-filter`: Regex for descriptions (capture group 1 appended)
- `--no-sort`: Disable sorting by name
- `--mincase` (default: 1): Min donors required
- `--case-col` (default: "case"): Donor column

**Example:**
```bash
immunediscover table fasta blast_V.tsv.gz novel_V.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 2 --cleanup " Novel"
```

---

### table collect

**Purpose:** Concatenate multiple TSVs (schema must match).

**Synopsis:**
```bash
immunediscover table collect <pattern> <output>
```

**Arguments:**
- `pattern`: Glob pattern (e.g., "plate*_exact.tsv.gz")
- `output`: Output file

**Example:**
```bash
immunediscover table collect "plate*_exact.tsv.gz" all_exact.tsv.gz
```

---

### table exclude

**Purpose:** Exclude rows matching reference FASTA by name or sequence.

**Synopsis:**
```bash
immunediscover table exclude <input> <output> <fasta> [options]
```

**Arguments:**
- `-n, --colname` (default: "allele_name"): Name column
- `-s, --colseq` (default: "seq"): Sequence column

---

## fasta merge

**Purpose:** Merge multiple FASTA files, keeping unique sequences.

**Synopsis:**
```bash
immunediscover fasta merge <output> <input1> <input2> [...] [options]
```

**Arguments:**
- `output`: Output merged FASTA
- `inputs`: Input FASTA files (2+)
- `-c, --cleanup`: Regex to remove from names
- `--no-sort`: Disable sorting
- `--prefer-last`: Prefer last name for duplicates (default: first)
- `--add-source-prefix`: Add filename prefix to names

**Example:**
```bash
immunediscover fasta merge combined.fasta known.fasta novel.fasta --cleanup " Novel"
```

---

## fasta diff

**Purpose:** Compare FASTA files by sequence identity.

**Synopsis:**
```bash
immunediscover fasta diff <fasta1> <fasta2> [...]
```

Prints union, intersection, and differences to stdout.

---

## fasta hash

**Purpose:** Add hash-based _S suffix to allele names.

**Synopsis:**
```bash
immunediscover fasta hash <input> > <output>
```

Outputs to stdout with names: `{name}_S{hash}`

