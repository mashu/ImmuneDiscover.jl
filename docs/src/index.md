```@meta
CurrentModule = immunediscover
```

[![Release](https://gitlab.com/gkhlab/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/gkhlab/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gkhlab.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/gkhlab/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/gkhlab/immunediscover.jl/commits/main)

# Immunediscover

Immunediscover is a comprehensive software package for analyzing genomic next-generation sequencing (NGS) data from immune receptors, specifically designed for discovering and characterizing immunoglobulin alleles. The package provides a complete pipeline from demultiplexing raw FASTQ reads to discovering novel alleles and performing association analysis.

## Overview

Immunediscover comprises multiple command groups, each containing specialized tools for different stages of analysis:

- **demultiplex** (top-level): Index-based demultiplexing of plate libraries
- **search**: Gene discovery tools (exact, blast, hsmm, heptamer)
- **analyze**: Quality control and analysis (association, haplotype, bwa)
- **table**: Data manipulation utilities
- **fasta**: FASTA file operations

The typical workflow starts with a FASTQ file and an index file, followed by demultiplexing, exact matching or discovery, and optional downstream analysis.

# Installation

## Downloading the Software

Immunediscover is provided as a pre-compiled binary for Linux systems. Download the appropriate version from the [releases page](https://gitlab.com/gkhlab/immunediscover.jl/-/releases).

## Setup Instructions

1. **Extract the Downloaded File**: Unzip the archive to your preferred directory (e.g., `/home/user/immunediscover`).

2. **Set Up Environment Variable**: Add the `bin` directory to your `PATH`:
   ```bash
   export PATH=$PATH:/home/user/immunediscover/bin/
   ```

3. **Make Change Permanent**: Add the export command to your `~/.bashrc` file:
   ```bash
   echo 'export PATH=$PATH:/home/user/immunediscover/bin/' >> ~/.bashrc
   ```

## Multithreading

Many subcommands use multithreading for improved performance. Control the number of threads with the `JULIA_NUM_THREADS` environment variable:

```bash
export JULIA_NUM_THREADS=8
immunediscover search blast input.tsv.gz database.fasta output.tsv.gz -g V
```

# Basic Usage Workflow

A typical analysis pipeline consists of:

1. **Demultiplexing**: Segregate reads by well/donor using plate indices
2. **Search**: Find alleles using exact matching, BLAST, or HSMM
3. **Analysis**: Perform association analysis, haplotype inference, or BWA QC
4. **Export**: Convert results to FASTA for downstream use

---

# Command Reference

## demultiplex

Demultiplex indexed plate libraries into per-read metadata.

### Purpose
Segregates FASTQ reads by forward and reverse barcodes, assigns to cases (donors), and filters by minimum read length.

### Usage
```bash
immunediscover demultiplex input.fastq.gz indices.tsv output.tsv.gz
```

### Input Files

**input.fastq.gz** (or `.fastq`): Single-end FASTQ file with sequencing reads.

**indices.tsv**: Tab-separated file with three required columns (no other columns allowed):
- `forward_index`: Forward barcode sequence (regex pattern)
- `reverse_index`: Reverse barcode sequence (regex pattern)
- `case`: Donor/case identifier (string)

### Parameters

**Positional:**
- `fastq` (required): Input FASTQ file with reads (gzipped or uncompressed)
- `indices` (required): TSV file with well indices/barcodes and case identifiers
- `output` (required): Output TSV (automatically compressed to .gz)

**Optional:**
- `-l, --length` (default: 200): Minimum read length to keep
- `-s, --split`: Write per-case FASTQ files to `{fastq}_split/` directory
- `-f, --forwardarrayindex` (default: ""): Name of forward array index for dual indexing
- `--case-filter-regex`: Regex pattern to filter cases (e.g., `"[ACDERF]"` keeps only matching cases)

### Output Files

**{output}.tsv.gz**: Demultiplexed reads with columns:
- `well`: Integer identifying the well (1-based index in indices file)
- `case`: Donor/case identifier
- `name`: Original FASTQ read identifier
- `genomic_sequence`: Full read sequence

**{output}.tsv.gz.log**: Per-well statistics with columns:
- `well`, `case`: Identifiers
- `len_μ`, `len_σ`: Mean and standard deviation of read lengths
- `len_q1`, `len_q2`, `len_q3`: Length quartiles
- `len_min`, `len_max`: Length range
- `seq_count`: Number of reads

### Example
```bash
# Basic demultiplexing
immunediscover demultiplex plate1.fastq.gz plate1_indices.tsv plate1_demux.tsv.gz

# With length filter and case filtering
immunediscover demultiplex plate1.fastq.gz plate1_indices.tsv plate1_demux.tsv.gz \
  --length 250 --case-filter-regex "[ACDERF]"

# Batch processing multiple plates
for fastq in *.fastq.gz; do
  base=${fastq%.fastq.gz}
  immunediscover demultiplex $fastq ${base}_indices.tsv ${base}_demux.tsv.gz
done
```

---

## search exact

Exact match search against a database of known alleles with robust filtering.

### Purpose
Identifies reads with exact matches to reference alleles, extracts flanking RSS (Recombination Signal Sequences) elements, and applies count and frequency filters to remove sequencing errors.

### Usage
```bash
immunediscover search exact demux.tsv.gz database.fasta output.tsv.gz -g V -c 5 -f 0.1
```

### Input Files

**demux.tsv.gz**: Demultiplexed TSV with columns `well`, `case`, `name`, `genomic_sequence`.

**database.fasta**: Reference allele database in FASTA format. Sequence names should follow format `GENE*ALLELE` (e.g., `IGHV1-69*01`).

### Parameters

**Positional:**
- `tsv` (required): TSV file with demultiplexed data
- `fasta` (required): FASTA file with query alleles
- `output` (required): TSV file to save output

**Required:**
- `-g, --gene`: Gene type - must be "V", "D", or "J"

**Optional:**
- `-c, --mincount` (default: 5): Minimum read count for an allele
- `-f, --minratio` (default: 0.1): Minimum allelic ratio within each gene group (minor/major allele)
- `--min-allele-mratio` (default: 0.05): Minimum ratio of allele count to cross-case median count
- `--min-gene-mratio` (default: 0.05): Minimum ratio of gene count to cross-case median count
- `--rss` (default: "heptamer"): Comma-separated RSS elements to extract: heptamer, spacer, nonamer
- `--extension`: Length of extension on RSS side (replaces --rss option)
- `-a, --affix` (default: 13): Bases to extract from non-RSS side
- `-t, --top` (default: 1): Maximum flank variants per allele/case (1 = collapsed mode)
- `-r, --refgene`: Reference gene(s) for ratio computation (space-separated)
- `-l, --limit` (default: 0): Limit input reads (0 = no limit, useful for testing)
- `-e, --expect`: TSV file with per-gene/allele ratio thresholds (columns: name, ratio)
- `-d, --deletion`: TSV file with per-gene deletion frequency thresholds (columns: name, ratio)
- `--locus` (default: "IG"): Locus prefix to filter genes for frequency calculations
- `--ref-fasta`: Optional reference FASTA to mark known vs novel sequences (adds `isin_db` column)
- `--raw`: Path to save unfiltered raw results for diagnostics
- `-n, --noplot`: Disable unicode boxplot of gene counts

### Output Files

**{output}.tsv.gz**: Exact matches with columns:
- `well`: Well identifier
- `case`: Donor/case identifier
- `gene`: Gene name (extracted from db_name by splitting on "*")
- `db_name`: Reference allele name
- `count`: Number of reads matching this allele sequence (collapsed)
- `full_count`: Number of reads matching sequence + flanks (uncollapsed)
- `ratio`: Allelic ratio within gene (count / max count in gene)
- `full_ratio`: Allelic ratio for full records
- `sequence`: Core allele sequence
- `prefix`: N bases before gene (affix)
- `suffix`: N bases after gene (affix) OR extension sequence
- `heptamer`, `spacer`, `nonamer`: RSS elements (if --rss specified)
- `flank_index`: Flank variant number (1 to top)
- `allele_freq`: Allele frequency within gene (per case)
- `allele_case_freq`: Allele frequency within case (all genes)
- `gene_case_freq`: Gene frequency within case
- `gene_count`: Total reads for this gene in case
- `case_count`: Total reads in case
- `cross_case_median_count`: Median count for this allele across all cases
- `cross_case_median_allele_count`: Median allele count across cases
- `cross_case_median_gene_count`: Median gene count across cases
- `allele_to_cross_case_median_ratio`: Allele count / cross-case median
- `gene_to_cross_case_median_ratio`: Gene count / cross-case median
- `isin_db`: "Novel" or "" (if --ref-fasta provided)
- `count_{refgene}_ratio`: Ratio to reference gene count (if --refgene specified)
- `gene_count_{refgene}_ratio`: Gene count ratio to reference (if --refgene specified)

### Examples

```bash
# Basic V gene search
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V

# D gene with extended flanks
immunediscover search exact demux.tsv.gz IGHD.fasta exact_D.tsv.gz -g D \
  --extension 40 --rss ""

# Low-stringency search for rare alleles
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V_lowstring.tsv.gz \
  -g V -c 1 -f 0.01 --min-allele-mratio 0.01

# With reference gene ratios
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V \
  --refgene "IGHV1-69*01" "IGHV3-23*01"

# Mark novel sequences
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V \
  --ref-fasta IGHV_reference.fasta
```

### Notes
- **Filtering logic**: Full records (sequence + flanks) are filtered first, then collapsed records
- **RSS extraction**: For V genes, RSS is 3' of gene; for J genes, RSS is 5' of gene; for D genes, RSS on both sides
- **Cross-case filtering**: Removes alleles with abnormally low counts compared to other cases (likely errors)
- **Expect/deletion files**: Override default ratio thresholds for specific genes (useful for known pseudogenes or deletions)

---

## search blast

BLAST-based allele discovery with sequence extension, trimming, and identity clustering.

### Purpose
Discovers novel alleles using BLAST alignment, with gene-specific presets for V, D, and J genes. Performs sequence extension (for D genes), re-alignment to database, and filtering based on coverage, frequency, and distance.

### Usage
```bash
immunediscover search blast demux.tsv.gz database.fasta output.tsv.gz -g D
```

### Input Files

**demux.tsv.gz**: Demultiplexed TSV with columns `well`, `case`, `name`, `genomic_sequence`.

**database.fasta**: Reference allele database.

### Parameters

**Positional:**
- `input` (required): TSV file with demultiplex data
- `fasta` (required): FASTA file with database sequences
- `output` (required): TSV file to save discovery results

**Gene Selection:**
- `-g, --gene`: Use gene preset (V, D, or J) - automatically sets optimal parameters
- `-G, --show-presets`: Display preset parameters for V, D, J and exit

**Extension (for D genes primarily):**
- `--forward` (default: 20): Forward extension length in nucleotides
- `--reverse` (default: 20): Reverse extension length in nucleotides

**BLAST Settings:**
- `-a, --args` (default: "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"): Additional blastn arguments
- `-d, --maxdist` (default: 20): Maximum mismatches allowed for alleles
- `-e, --edge` (default: 0): Minimum nucleotides required between gene and read end
- `-s, --subjectcov` (default: 0.1): Minimum subject (database) coverage fraction

**Filtering:**
- `-c, --minfullcount` (default: 5): Minimum full cluster size (before collapsing)
- `-f, --minfullfreq` (default: 0.1): Minimum allelic frequency within gene group
- `-l, --length` (default: 290): Minimum length of aligned read
- `-q, --minquality` (default: 0.75): Minimum alignment quality for trimming (mismatch/length)
- `--min-corecov` (default: 0.6): Minimum ratio of aligned length to database sequence length
- `-i, --isin`: Keep truncated substrings of known alleles (default: filter them out)
- `--keep-failed`: Keep rows where trimming failed (default: drop them)

**Other:**
- `-p, --pseudo`: FASTA file with pseudo-genes (prefixed with "P" in combined database)
- `-o, --overwrite`: Overwrite existing BLAST cache files
- `-v, --verbose`: Save intermediate files and print detailed output

### Gene Presets

**V Gene Preset:**
- forward: 12, reverse: 12
- minfullfreq: 0.1
- length: 290
- maxdist: 10
- minfullcount: 10
- args: "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"

**D Gene Preset:**
- forward: 40, reverse: 40
- minfullfreq: 0.2
- length: 5
- maxdist: 20
- minfullcount: 10
- edge: 10
- subjectcov: 0.25
- minquality: 0.5
- args: "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5"

**J Gene Preset:**
- forward: 12, reverse: 12
- minfullfreq: 0.1
- length: 10
- maxdist: 10
- minfullcount: 10
- args: "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"

### Output Files

**{output}.tsv.gz**: Discovered alleles with columns:
- `well`, `case`: Identifiers
- `qseqid`: Query sequence ID
- `sseqid`: Subject (database) sequence ID
- `gene`: Gene name (from sseqid)
- `qseq`: Original query sequence
- `db_seq`: Database sequence
- `prefix`, `suffix`: Extension sequences
- `aln_qseq`: Aligned query sequence (gaps removed, trimmed)
- `aln_mismatch`: Mismatches in aligned core
- `corecov`: Coverage of database gene (length(aln_qseq)/length(db_seq))
- `isin_db`: Whether aligned sequence is substring of known allele
- `full_count`: Cluster size
- `full_frequency`: Frequency within well/case/gene
- `count`: Aggregated count for this allele+flanks
- `frequency`: Frequency after aggregation
- `allele_name`: Assigned name (original or with "_S{hash} Novel" suffix)
- BLAST output columns: `pident`, `nident`, `length`, `mismatch`, `gapopen`, `qcovs`, `qcovhsp`, `qstart`, `qend`, `sstart`, `send`, `qlen`, `slen`, `evalue`, `bitscore`, `sstrand`

**{fasta}-combined.fasta**: Combined database (pseudo-genes + regular genes)

**{fasta}-combined-extended.fasta**: Extended sequences (if forward/reverse > 0)

**{fasta}.affixes**: TSV with extension information (columns: name, prefix, suffix)

**{input}-blast.tsv**: BLAST raw output (if verbose)

### Workflow

1. **Extension** (optional): Extends database sequences by most common flanks from reads
2. **BLAST**: Aligns reads to extended database
3. **Re-alignment**: Trims extensions and re-aligns to original database using Needleman-Wunsch
4. **Filtering**: Applies coverage, quality, frequency, and distance filters
5. **Naming**: Assigns allele names with hash suffix for novel sequences

### Examples

```bash
# V gene discovery with preset
immunediscover search blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz -g V

# D gene discovery with preset
immunediscover search blast demux_D.tsv.gz IGHD.fasta blast_D.tsv.gz -g D

# Custom parameters overriding preset
immunediscover search blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz \
  -g V --forward 20 --maxdist 5

# With pseudo-genes
immunediscover search blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz \
  -g V --pseudo IGHV_pseudo.fasta

# Low-stringency for rare alleles
immunediscover search blast demux.tsv.gz IGHV.fasta blast_lowstring.tsv.gz \
  -g V -c 1 -f 0.05 --min-corecov 0.5

# View presets
immunediscover search blast --show-presets
```

### Notes
- **Extension purpose**: For D genes, extension captures flanking V/J sequences to improve BLAST sensitivity
- **Cache**: BLAST results are cached; use `--overwrite` to regenerate
- **Quality threshold**: Alignment quality is `(length - mismatches) / length`
- **Novel allele naming**: Format is `{closest_gene}_S{hash} Novel` where hash is MD5-based

---

## search hsmm

Detect D genes using a Hidden Semi-Markov Model trained on RSS flanks.

### Purpose
Discovers D genes by training an HSMM on RSS (Recombination Signal Sequences) flanks from known D alleles, then applying the model to all reads to detect D segments with both known and novel sequences. 

### Usage
```bash
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz
```

### Input Files

**demux.tsv.gz**: Demultiplexed TSV with columns `well`, `case`, `name`, `genomic_sequence`.

**IGHD.fasta**: Known D allele database for training.

### Parameters

**Positional:**
- `tsv` (required): TSV/TSV.GZ demultiplex file
- `fasta` (required): FASTA file with D alleles (known reference)
- `output` (required): TSV.GZ file to save detected D alleles

**Training Filters:**
- `-r, --ratio` (default: 0.2): Allelic ratio threshold for known D selection per donor and gene
- `-c, --mincount` (default: 10): Minimum count for a known D allele to be considered in training

**Model Parameters:**
- `--min-gene-len` (default: 10): Minimum D gene length for HSMM duration model (auto if 0)
- `--max-gene-len` (default: 70): Maximum D gene length for HSMM duration model (auto if 0)

**Detection Filters:**
- `--min-posterior` (default: 0.7): Minimum posterior probability for accepting HSMM detection
- `--min-heptamer-prob-pre` (default: 0.05): Minimum probability under pre-heptamer PWM (0 disables)
- `--min-heptamer-prob-post` (default: 0.05): Minimum probability under post-heptamer PWM (0 disables)

**Output Filters:**
- `--out-mincount` (default: 10): Minimum count for extracted D to keep in output
- `--out-minratio` (default: 0.2): Minimum allelic ratio within gene per donor for output

**Other:**
- `-l, --limit` (default: 0): Limit demultiplexed reads to process (0 = no limit)

### Output Files

**{output}.tsv.gz**: Detected D genes with columns:
- `well`, `case`: Identifiers
- `sequence`: Detected D gene sequence
- `pre_nonamer`, `pre_spacer`, `pre_heptamer`: 5' RSS elements (9+12+7 bp)
- `post_heptamer`, `post_spacer`, `post_nonamer`: 3' RSS elements (7+12+9 bp)
- `heptamer_logp_pre`, `heptamer_logp_post`: Log-probabilities under heptamer PWMs
- `heptamer_prob_pre`, `heptamer_prob_post`: Probabilities (exp of log-probabilities)
- `log_path_prob`: Log-probability of best path
- `log_total_prob`: Log-probability summed over all paths
- `posterior_prob`: Posterior probability of best path (exp(log_path_prob - log_total_prob))
- `isin_db`: Boolean indicating if sequence is in reference database
- `db_name`: Reference allele name if exact match, "Novel" otherwise
- `nearest_db`: Name of closest reference allele (by Levenshtein distance)
- `nearest_db_dist`: Levenshtein distance to nearest reference
- `allele_name`: Final allele name (reference name or with _S{hash} suffix)
- `count`: Number of reads with this D sequence (per well/case)
- `ratio`: Allelic ratio within gene group (per well/case)
- `gene`: Gene name (from nearest_db)

### Algorithm

1. **Known D search**: Exact search for known D alleles in reads, extract RSS flanks
2. **Training set selection**: Per donor/gene, keep top 2 alleles by count and ratio ≥ threshold
3. **HSMM training**: Build position weight matrices (PWMs) for RSS elements and IID model for D gene from training flanks
4. **Detection**: Apply HSMM to all reads, find best-scoring D segment with RSS flanks
5. **Filtering**: Filter by posterior probability, heptamer PWM probabilities, count, and ratio

### HSMM Structure

The model consists of 7 sequential states:
- pre_nonamer (9 bp, PWM)
- pre_spacer (12 bp, PWM)
- pre_heptamer (7 bp, PWM)
- **D gene (variable length, IID emission + duration model)**
- post_heptamer (7 bp, PWM)
- post_spacer (12 bp, PWM)
- post_nonamer (9 bp, PWM)

### Examples

```bash
# Basic D gene detection
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz

# Relaxed training filter for rare alleles
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --ratio 0.1 --mincount 5

# Strict detection thresholds
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --min-posterior 0.9 --min-heptamer-prob-pre 0.1 --min-heptamer-prob-post 0.1

# Custom gene length range
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz \
  --min-gene-len 8 --max-gene-len 50

# Test on subset
immunediscover search hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz --limit 10000
```

### Notes
- **Training quality**: Requires high-quality known D alleles across multiple donors
- **Posterior probability**: Ratio of best path probability to sum of all paths (confidence metric)
- **Heptamer filters**: PWM-based quality filters for RSS elements (0 disables)
- **Parallelization**: Detection phase is parallelized using multiple threads
- **Levenshtein distance**: Used for assigning closest known allele to novel sequences

---

## search heptamer

Identify heptamer RSS positions and extend/trim V reads accordingly.

### Purpose
Locates heptamer recombination signal sequences (RSS) downstream of V genes, extends read sequences to include them, and summarizes results by collapsing identical sequences.

### Usage
```bash
immunediscover search heptamer demux.tsv.gz IGHV.fasta heptamer.tsv.gz summary.tsv --chain IGHV
```

### Input Files

**demux.tsv.gz**: Demultiplexed TSV with columns `well`, `case`, `name`, `genomic_sequence`.

**IGHV.fasta**: V gene reference database.

**heptamers.json**: JSON file with heptamer sequences per chain (auto-generated if missing).

### Parameters

**Positional:**
- `tsv` (required): TSV file with demultiplexed reads
- `fasta` (required): FASTA file with query alleles
- `output` (required): TSV file to save data with identified heptamers
- `summary` (required): TSV file to save summary collapsed alleles with statistics

**Optional:**
- `-j, --json` (default: "heptamers.json"): JSON file with heptamer dictionary
- `-c, --chain` (default: "IGHV"): Chain - must be IGKV, IGLV, or IGHV
- `-d, --maxdist` (default: 1): Maximum Hamming distance from heptamers in JSON
- `-b, --begin` (default: 0): Bases to trim from 5' beginning of query sequence
- `-e, --end` (default: 8): Bases to trim from 3' end of query sequence
- `-m, --mincount` (default: 1): Minimum count allowed in summary
- `-r, --ratio` (default: 0.25): Minimum ratio between full allele count and trimmed allele count

### Heptamer Sequences

If `heptamers.json` doesn't exist, default heptamers are created:
- **IGHV**: CACAGTG, CACAATG, CACAGAG, CACGGTG, CACAGCG
- **IGKV**: CACAGTG, CACACTG, CACTGTG, CACGGTG, CACAATG, CACATTG
- **IGLV**: CACAGTG, CACGGTG, CATGGTG, CACGCTG, CACAGCG, CACAGTA, CATAGTG, CACAATG

### Output Files

**{output}.tsv.gz**: Per-read heptamer identifications with columns:
- `well`, `case`, `name`, `genomic_sequence`: From input
- `db_name`: Reference allele name
- `full_match`: Range of full database sequence in read (if found)
- `trimmed_match`: Range of trimmed database sequence in read
- `db_length`: Length of database sequence
- `full_length`: Length of extended sequence (trimmed + heptamer extension)
- `full_sequence`: Extended sequence
- `heptamer`: Identified heptamer sequence (7 bp)

**{summary}.tsv**: Collapsed summary with columns:
- `db_name`: Original reference allele name
- `new_name`: Assigned name for novel alleles (_S{hash} suffix)
- `full_count`: Number of reads with exact match to full (untrimmed) database sequence
- `trimmed_count`: Number of reads with exact match to trimmed database sequence
- `case`: Comma-separated list of cases with this allele
- `allele_count`: Number of reads with this full_sequence + heptamer combination
- `ratio`: allele_count / trimmed_count
- `full_length`, `db_length`: Sequence lengths
- `sequence`: Full extended sequence
- `heptamer`: Heptamer sequence
- `suffix`: Sequence downstream of trimmed match (includes heptamer)

### Examples

```bash
# Basic heptamer search for IGHV
immunediscover search heptamer demux.tsv.gz IGHV.fasta heptamer.tsv.gz summary.tsv \
  --chain IGHV

# Allow 2 mismatches in heptamer
immunediscover search heptamer demux.tsv.gz IGHV.fasta heptamer.tsv.gz summary.tsv \
  --chain IGHV --maxdist 2

# Trim 5 bp from start, 10 bp from end
immunediscover search heptamer demux.tsv.gz IGHV.fasta heptamer.tsv.gz summary.tsv \
  --chain IGHV --begin 5 --end 10

# Strict filtering in summary
immunediscover search heptamer demux.tsv.gz IGHV.fasta heptamer.tsv.gz summary.tsv \
  --chain IGHV --mincount 10 --ratio 0.5
```

### Notes
- **Trimming**: begin/end parameters trim query before searching (removes primer regions)
- **Heptamer search**: Prioritizes lower Hamming distances, then earlier positions in suffix
- **Ratio filter**: Ensures extended sequence is common enough relative to trimmed match

---

## analyze association

Analyze cross-donor association (co-occurrence) among alleles.

### Purpose
Computes pairwise association metrics (phi coefficient) between alleles across multiple donors/cases to identify alleles that co-occur (potential haplotype pairs) or anti-correlate. Performs hierarchical clustering to group associated alleles.

### Methodology Note

This analysis computes **phi coefficient** (correlation-based association) rather than **classical linkage disequilibrium (LD)** measures because:

**Our data is unphased**: We know which alleles are present in each donor, but not which chromosome (maternal vs. paternal) they're on. Classical LD measures like D (= pAB - pA×pB) and D' require phased haplotype data showing which specific alleles are physically together on the same chromosome.

**What we compute:**
- **Phi coefficient (φ)**: Pearson's correlation for binary presence/absence data across donors
- **φ²** (r²): Mathematically equivalent to the standardized LD r² used in population genetics
- Measures: "Do alleles A and B co-occur across donors more than expected by chance?"

**Equivalence**: Using `--similarity r2` provides results comparable to classical LD r² analysis, but computed from unphased genotype data rather than requiring phased haplotypes.

### Usage
```bash
immunediscover analyze association input.tsv.gz edges.tsv.gz --case-col case --allele-col db_name
```

### Input Files

**input.tsv.gz**: Any TSV with at least two columns (names configurable):
- Case/donor identifier column (default: `case`)
- Allele identifier column (default: `db_name`)

### Parameters

**Positional:**
- `input` (required): TSV/TSV.GZ with columns for case and allele
- `edges` (required): TSV.GZ to save pairwise edges table

**Column Names:**
- `-C, --case-col` (default: "case"): Name of donor/case id column
- `-A, --allele-col` (default: "db_name"): Name of allele column

**Filtering:**
- `-m, --min-donors` (default: 2): Minimum donors required to include an allele
- `--min-support` (default: 3): Minimum co-present donors (n11) to consider an edge
- `--min-jaccard` (default: 0.2): Minimum Jaccard co-presence to consider an edge

**Clustering:**
- `--similarity` (default: "r"): Similarity mode - "r" (positive correlation) or "r2" (LD-style)
- `--threshold` (default: 0.5): Similarity threshold for complete-linkage clustering
- `--min-cluster-size` (default: 3): Minimum cluster size to output
- `--clusters`: Optional path to save clusters TSV

### Output Files

**{edges}.tsv.gz**: Pairwise association metrics with columns:
- `allele_a`, `allele_b`: Pair of alleles
- `r`: Phi coefficient (Pearson correlation for binary presence/absence data)
- `r2`: r squared (proportion of variance explained)
- `jaccard`: Jaccard index = n11 / (n11 + n10 + n01)
- `support`: n11 (number of donors with both alleles)
- `similarity`: Computed similarity (r or r2, filtered by jaccard and support)

**{clusters}.tsv** (if --clusters specified): Cluster assignments with columns:
- `group_id`: Cluster number
- `allele`: Allele name
- `donors`: Comma-separated list of donors with this allele
- `n_donors`: Number of donors
- `mean_n11`, `mean_n10`, `mean_n01`, `mean_n00`: Average pairwise contingency counts within cluster
- `mean_r`, `max_r`, `min_r`: Association statistics within cluster (phi coefficient values)

### Association Metrics

**Contingency Table** (for alleles A and B):
- `n11`: Donors with both A and B
- `n10`: Donors with A but not B
- `n01`: Donors with B but not A
- `n00`: Donors with neither A nor B

**Phi Coefficient (r)** - Pearson correlation for binary variables:
```
r = (n11*n00 - n10*n01) / sqrt((n11+n10)(n11+n01)(n10+n00)(n01+n00))
```
This is equivalent to Pearson's r for 2×2 contingency tables, measuring correlation between allele presence/absence patterns across donors. Ranges from -1 (perfect negative association) to +1 (perfect positive association).

**Relationship to classical LD:**
- φ² (phi squared, r²) is mathematically equivalent to the standardized r² used in LD analysis
- However, classical LD coefficient D = pAB - pA×pB requires phased haplotypes
- Since our data is unphased (we don't know which alleles are on the same chromosome), we use the correlation-based phi coefficient
- The `--similarity r2` option provides an LD-comparable metric without requiring phase information

**Jaccard Index**:
```
J = n11 / (n11 + n10 + n01)
```
Measures overlap ignoring double-absences (n00), ranging from 0 to 1.

### Examples

```bash
# Basic association analysis
immunediscover analyze association exact_V.tsv.gz association_edges.tsv.gz

# With custom column names
immunediscover analyze association input.tsv.gz edges.tsv.gz \
  --case-col donor_id --allele-col allele_name

# Strict filtering for high-confidence associations
immunediscover analyze association exact_V.tsv.gz edges.tsv.gz \
  --min-donors 5 --min-support 10 --min-jaccard 0.5

# r² similarity with clustering output
immunediscover analyze association exact_V.tsv.gz edges.tsv.gz \
  --similarity r2 --threshold 0.7 --clusters association_clusters.tsv

# Low stringency for exploratory analysis
immunediscover analyze association exact_V.tsv.gz edges.tsv.gz \
  --min-donors 1 --min-support 1 --min-jaccard 0.1
```

### Notes
- **r vs r²**: "r" mode keeps only positive correlations (haplotype pairs); "r²" includes negative correlations, and is equivalent to standardized LD r²
- **Unphased data**: This analysis measures population-level co-occurrence without requiring phase information about which alleles are on the same chromosome
- **Phased vs unphased**: If you had phased haplotype data, classical LD measures (D, D') would be more appropriate; for unphased genotype data (which alleles are present per donor), phi coefficient is the standard choice
- **Filtering**: Edges filtered by min-support and min-jaccard before similarity computation
- **Clustering**: Complete-linkage hierarchical clustering on distance = 1 - similarity
- **Use cases**: Identifying haplotype blocks, validating novel alleles, detecting sequencing artifacts

---

## analyze haplotype

Infer approximate haplotypes per case using diploid assumptions.

### Purpose
Infers homozygous vs heterozygous genotypes for each donor and gene based on allelic ratios. Assumes diploid genetics (≤2 functional alleles per gene per donor).

### Usage
```bash
immunediscover analyze haplotype input.tsv.gz haplotypes.tsv
```

### Input Files

**input.tsv.gz**: TSV with columns (names configurable):
- Case/donor identifier (default: `case`)
- Allele identifier (default: `db_name`)
- Gene identifier (default: `gene`)
- `count`: Read count per allele (required)

Typically uses output from `search exact` or `search blast`.

### Parameters

**Positional:**
- `input` (required): TSV/TSV.GZ with allele data
- `output` (required): TSV file to save haplotype results

**Column Names:**
- `-C, --case-col` (default: "case"): Donor/case identifier column
- `-A, --allele-col` (default: "db_name"): Allele name column
- `-G, --gene-col` (default: "gene"): Gene name column

**Filtering:**
- `-c, --mincount` (default: 5): Minimum count for allele to be considered
- `-r, --min-ratio` (default: 0.1): Minimum minor/major allele ratio for heterozygous classification

**Optional:**
- `-f, --novel-fasta`: FASTA file with novel alleles (adds novel_1, novel_2 columns)

### Output Files

**{output}.tsv**: Haplotype inferences with columns:
- `case`: Donor identifier
- `gene`: Gene identifier
- `genotype`: "homozygous", "heterozygous", or "uncertain"
- `allele_1`: Primary allele (highest count)
- `allele_2`: Secondary allele (second highest count, empty if homozygous)
- `count_1`: Count for primary allele
- `count_2`: Count for secondary allele (0 if homozygous)
- `ratio`: Minor/major allele ratio (count_2 / count_1)
- `total_count`: Total count for this gene in this donor
- `other_alleles`: Comma-separated additional alleles (if >2 present)
- `novel_1`, `novel_2`: Boolean indicating if allele is novel (if --novel-fasta provided)

### Classification Logic

**Homozygous**: One allele present OR ratio < min-ratio
**Heterozygous**: Two alleles with ratio ≥ min-ratio
**Uncertain**: More than two alleles present (unusual for diploid)

### Examples

```bash
# Basic haplotype inference
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv

# Mark novel alleles
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv \
  --novel-fasta novel_V.fasta

# Relaxed heterozygosity threshold
immunediscover analyze haplotype exact_V.tsv.gz haplotypes.tsv \
  --min-ratio 0.05 --mincount 3

# Custom column names
immunediscover analyze haplotype input.tsv.gz haplotypes.tsv \
  --case-col donor --allele-col allele_name --gene-col gene_id
```

### Notes
- **Diploid assumption**: May not hold for pseudogenes or copy number variants
- **Ratio interpretation**: 0.5 = balanced heterozygous; <0.1 likely homozygous or artifact
- **Novel alleles**: Useful for assessing how many novel alleles are present in heterozygous genotypes

---

## analyze bwa

Genome mapping QC using BWA to retain sequences mapped to target chromosome.

### Purpose
Aligns candidate allele sequences to a reference genome using BWA-MEM, filters for target chromosome, and adds alignment metadata. Used for quality control to remove off-target sequences and pseudogenes.

### Usage
```bash
immunediscover analyze bwa input.tsv.gz output.tsv.gz genome.fasta --chromosome "chromosome 14"
```

### Input Files

**input.tsv.gz**: TSV with columns:
- Allele name column (default: `best_name`)
- Sequence column(s) (default: `prefix`, `best_aln`, `suffix` - concatenated)

**genome.fasta** (or multiple): BWA-indexed reference genome. Index with:
```bash
bwa index genome.fasta
```

### Parameters

**Positional:**
- `tsv` (required): TSV file with allele_name and seq columns
- `genome` (required, 1+): FASTA file(s) with indexed genome
- `output` (required): TSV file to save filtered input

**Optional:**
- `-c, --chromosome` (default: "chromosome 14"): Chromosome string to filter by
- `-t, --tag` (default: "(.*Primary Assembly.*)|(.*alternate locus.*)"): Regex to filter valid chromosome descriptions
- `-n, --colname` (default: "best_name"): Column with allele names
- `-s, --colseq` (default: ["prefix", "best_aln", "suffix"]): Columns to concatenate for sequence

### Output Files

**{output}.tsv.gz**: Filtered rows (only sequences mapping to target chromosome) with added columns:
- `position`: Chromosome and position (format: "chr_name:pos")
- `edit_distance`: Number of mismatches in alignment (NM field)
- `ref_sequence`: Reference sequence at alignment position (orientation-adjusted)
- `orientation`: "+" for forward strand, "-" for reverse strand
- Original columns preserved

**/tmp/discarded.tsv**: List of allele names that mapped to off-target chromosomes

### Examples

```bash
# Basic BWA filtering for human chromosome 14
bwa index GCF_000001405.25.fasta
immunediscover analyze bwa candidates.tsv.gz filtered.tsv.gz \
  GCF_000001405.25.fasta --chromosome "chromosome 14"

# Mouse chromosome 12
immunediscover analyze bwa candidates.tsv.gz filtered.tsv.gz \
  mouse_genome.fasta --chromosome "chromosome 12"

# Custom column names
immunediscover analyze bwa candidates.tsv.gz filtered.tsv.gz \
  genome.fasta --colname allele_id --colseq sequence

# Concatenate multiple columns
immunediscover analyze bwa candidates.tsv.gz filtered.tsv.gz \
  genome.fasta --colseq prefix core suffix
```

### Notes
- **Indexing required**: Genome must be BWA-indexed before use
- **Multiple genomes**: Can provide multiple FASTA files (e.g., separate chromosomes)
- **Orientation**: ref_sequence is reverse-complemented if alignment is on reverse strand
- **Use case**: Essential QC step before validating novel alleles

---

## table Commands

Utilities for manipulating TSV files (join, filter, sort, aggregate, etc.).

### table outerjoin

Outer join two TSV files on specified key columns.

**Usage:**
```bash
immunediscover table outerjoin left.tsv.gz right.tsv.gz output.tsv.gz \
  --keys well,case
```

**Parameters:**
- `left`, `right`, `output` (required): File paths
- `-k, --keys` (required): Comma-separated key columns
- `--left-keys`, `--right-keys`: Different keys for each file
- `--left-prefix`, `--right-prefix`: Prefix for non-key columns
- `--left-select`, `--right-select`: Subset of columns to keep

**Output:** Combined TSV with all rows from both files, missing values filled with `missing`.

---

### table leftjoin

Left join two TSV files (keep all left rows).

**Usage:**
```bash
immunediscover table leftjoin left.tsv.gz right.tsv.gz output.tsv.gz \
  --keys case,gene
```

**Parameters:** Same as `outerjoin`.

**Output:** Combined TSV preserving all left file rows.

---

### table transform

Regex-based column transformation with capture groups.

**Usage:**
```bash
immunediscover table transform input.tsv.gz output.tsv.gz \
  --column db_name --pattern 'IGHV(.*)\\*(.*)' --replacement 'V\\1-\\2'
```

**Parameters:**
- `-c, --column` (required): Target column
- `-p, --pattern` (required): Regex with capture groups
- `-r, --replacement` (required): Replacement using \\1, \\2, etc.
- `--new-column`: Optional new column name (preserves original)

**Output:** TSV with transformed column.

---

### table aggregate

Group by columns and count unique groups.

**Usage:**
```bash
immunediscover table aggregate input.tsv.gz output.tsv.gz \
  --group-by case,gene,allele_name
```

**Parameters:**
- `-g, --group-by` (required): Comma-separated grouping columns
- `-k, --keep-columns`: Additional columns to keep (first value per group)
- `-c, --count-column` (default: "count"): Name for count column

**Output:** TSV with unique combinations and counts.

---

### table unique

Select distinct rows by specified columns.

**Usage:**
```bash
immunediscover table unique input.tsv.gz output.tsv.gz \
  --columns case,allele_name,sequence
```

**Parameters:**
- `-c, --columns` (required): Comma-separated columns for uniqueness

**Output:** TSV with duplicate rows removed.

---

### table sort

Sort TSV by one or more columns.

**Usage:**
```bash
immunediscover table sort input.tsv.gz output.tsv.gz \
  --columns case,count --reverse
```

**Parameters:**
- `-c, --columns` (required): Comma-separated sort columns (priority order)
- `-r, --reverse`: Sort descending (default ascending)

**Output:** Sorted TSV.

---

### table filter

Filter rows by regex pattern or numeric threshold.

**Usage:**
```bash
# Regex mode
immunediscover table filter input.tsv.gz output.tsv.gz \
  --column allele_name --pattern "Novel"

# Numeric mode
immunediscover table filter input.tsv.gz output.tsv.gz \
  --column count --operator ">=" --threshold 10
```

**Parameters:**
- `-c, --column` (required): Target column
- `--pattern`: Regex pattern (string filtering)
- `--operator`: Comparison operator (<, <=, >=, >)
- `--threshold`: Numeric value (used with --operator)

**Output:** Filtered TSV.

---

### table select

Project subset of columns.

**Usage:**
```bash
immunediscover table select input.tsv.gz output.tsv.gz \
  --columns case,gene,allele_name,count
```

**Parameters:**
- `-c, --columns` (required): Comma-separated columns to keep

**Output:** TSV with selected columns only.

---

### table fasta

Export sequences from TSV to FASTA.

**Usage:**
```bash
immunediscover table fasta input.tsv.gz output.fasta \
  --colname allele_name --colseq sequence --filter "Novel"
```

**Parameters:**
- `-n, --colname` (default: "allele_name"): Column with sequence names
- `-s, --colseq` (default: "seq"): Column with sequences
- `-d, --coldesc`: Optional column for FASTA descriptions
- `-f, --filter`: Regex to filter sequence names
- `-c, --cleanup`: Regex to remove from names (e.g., " Novel")
- `--desc-filter`: Regex for description (capture group 1 appended to header)
- `--no-sort`: Disable sorting by name
- `--mincase` (default: 1): Minimum donors required
- `--case-col` (default: "case"): Donor column name

**Output:** FASTA file with unique sequences.

---

### table collect

Concatenate multiple TSV files (schema must match).

**Usage:**
```bash
immunediscover table collect "plate*-exact.tsv.gz" collected.tsv.gz
```

**Parameters:**
- `pattern` (required): Glob pattern for input files
- `output` (required): Output file

**Output:** Combined TSV with all input rows.

---

### table exclude

Exclude rows if sequence name or sequence matches reference FASTA.

**Usage:**
```bash
immunediscover table exclude input.tsv.gz output.tsv.gz reference.fasta \
  --colname allele_name --colseq sequence
```

**Parameters:**
- `-n, --colname` (default: "allele_name"): Column with allele names
- `-s, --colseq` (default: "seq"): Column with sequences

**Output:** TSV with matching rows removed.

---

## fasta Commands

Operations on FASTA files.

### fasta merge

Merge multiple FASTA files, keeping only unique sequences.

**Usage:**
```bash
immunediscover fasta merge output.fasta input1.fasta input2.fasta input3.fasta
```

**Parameters:**
- `output` (required): Output merged FASTA
- `inputs` (required, 2+): Input FASTA files
- `-c, --cleanup`: Regex to remove from names
- `--no-sort`: Disable sorting by name
- `--prefer-last`: When duplicates have different names, prefer last (default: first)
- `--add-source-prefix`: Add filename prefix to sequence names

**Output:** FASTA with unique sequences. Duplicates reported to log.

---

### fasta diff

Compare FASTA files by sequence identity.

**Usage:**
```bash
immunediscover fasta diff file1.fasta file2.fasta file3.fasta
```

**Parameters:**
- `fasta` (required, 2+): FASTA files to compare

**Output:** Prints set operations (union, intersection, differences) to stdout.

---

### fasta hash

Add hash-based _S suffix to allele names.

**Usage:**
```bash
immunediscover fasta hash input.fasta > output.fasta
```

**Parameters:**
- `fastain` (required): Input FASTA

**Output:** FASTA with names in format `{name}_S{hash}` printed to stdout.

---

# Typical Workflows

## Workflow 1: V Gene Discovery

```bash
# 1. Demultiplex
immunediscover demultiplex plate_V.fastq.gz indices.tsv demux_V.tsv.gz --length 250

# 2. Exact search for known alleles
immunediscover search exact demux_V.tsv.gz IGHV_known.fasta exact_V.tsv.gz -g V

# 3. BLAST discovery of novel alleles
immunediscover search blast demux_V.tsv.gz IGHV_known.fasta blast_V.tsv.gz -g V

# 4. Export novel alleles to FASTA
immunediscover table fasta blast_V.tsv.gz novel_V.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 2

# 5. BWA QC to filter off-target
bwa index GCF_000001405.25.fasta
immunediscover table fasta blast_V.tsv.gz candidates.fasta --colname allele_name --colseq aln_qseq
# Create temporary TSV for BWA
immunediscover table select blast_V.tsv.gz bwa_input.tsv.gz --columns allele_name,aln_qseq
immunediscover analyze bwa bwa_input.tsv.gz bwa_filtered.tsv.gz GCF_000001405.25.fasta

# 6. Merge with known database
immunediscover fasta merge IGHV_updated.fasta IGHV_known.fasta novel_V.fasta --cleanup " Novel"

# 7. Re-run exact search with updated database
immunediscover search exact demux_V.tsv.gz IGHV_updated.fasta exact_V_final.tsv.gz -g V

# 8. Association analysis
immunediscover analyze association exact_V_final.tsv.gz association_edges.tsv.gz \
  --min-donors 3 --clusters association_clusters.tsv

# 9. Haplotype inference
immunediscover analyze haplotype exact_V_final.tsv.gz haplotypes.tsv \
  --novel-fasta novel_V.fasta
```

## Workflow 2: D Gene Discovery

```bash
# 1. Demultiplex
immunediscover demultiplex plate_D.fastq.gz indices.tsv demux_D.tsv.gz --length 200

# 2. Option A: BLAST with extension
immunediscover search blast demux_D.tsv.gz IGHD_known.fasta blast_D.tsv.gz -g D

# 2. Option B: HSMM detection (better for short D genes)
immunediscover search hsmm demux_D.tsv.gz IGHD_known.fasta hsmm_D.tsv.gz

# 3. Export and validate
immunediscover table fasta hsmm_D.tsv.gz novel_D.fasta \
  --colname allele_name --colseq sequence --filter "Novel" --mincase 3

# 4. Merge databases
immunediscover fasta merge IGHD_updated.fasta IGHD_known.fasta novel_D.fasta
```

## Workflow 3: Multi-Plate Batch Processing

```bash
# Demultiplex all plates
for fastq in plate*.fastq.gz; do
  base=${fastq%.fastq.gz}
  immunediscover demultiplex $fastq ${base}_indices.tsv ${base}_demux.tsv.gz
done

# Exact search on all plates
for demux in *_demux.tsv.gz; do
  base=${demux%_demux.tsv.gz}
  immunediscover search exact $demux IGHV.fasta ${base}_exact.tsv.gz -g V
done

# Collect results
immunediscover table collect "plate*_exact.tsv.gz" all_plates_exact.tsv.gz

# Aggregate by allele across all plates
immunediscover table aggregate all_plates_exact.tsv.gz summary.tsv.gz \
  --group-by db_name,sequence --keep-columns gene
```

---

# Tips and Best Practices

## Parameter Selection

**Minimum Count Thresholds:**
- High confidence: `--mincount 10` (default for most commands)
- Moderate: `--mincount 5`
- Rare alleles: `--mincount 1` (risk of false positives)

**Frequency Thresholds:**
- Homozygous-heavy: `--minratio 0.2`
- Balanced: `--minratio 0.1` (default)
- Low-expression alleles: `--minratio 0.05`

**Distance Thresholds:**
- V genes: `--maxdist 10` (default)
- D genes: `--maxdist 20` (default, D genes more variable)
- Very stringent: `--maxdist 5`

## Performance Optimization

**Threading:**
```bash
export JULIA_NUM_THREADS=16  # Use all available cores
```

**Memory:**
- Large datasets: Split by plate, process separately
- Limit input: `--limit 100000` for testing

## Quality Control Checklist

1. Check demultiplex statistics (wells with very low counts may indicate barcode issues)
2. Verify gene counts are balanced (plot from exact search)
3. Run BWA QC on novel alleles before validation
4. Check association analysis for unexpected patterns (may indicate artifacts)
5. Validate haplotypes match expected heterozygosity rates
6. Compare cross-case medians to identify donor-specific errors

## Common Issues

**Low demultiplexing rate:**
- Check barcode file format (no extra columns)
- Verify barcode orientation matches library prep
- Consider adapter trimming if indices are offset

**Too many novel alleles:**
- Increase `--mincount` and `--minfullfreq`
- Check `--min-corecov` (low values allow partial matches)
- Run BWA QC to filter pseudogenes
- Increase `--minquality` for better alignment filtering

**Missing expected alleles:**
- Decrease `--mincount` temporarily for exploration
- Check `--length` filter (may exclude short alleles)
- Verify locus filtering (`--locus`) includes target genes
- Check demux read length distribution

**High memory usage:**
- Process plates individually, then collect
- Use `--limit` for testing
- Reduce `--top` parameter (fewer flank variants)

**Unexpected associations:**
- Check for plate effects or batch artifacts
- Verify donors are truly independent
- Look for contamination or barcode swapping

---

# Output Column Glossary

## Common Columns (appear in multiple outputs)

- `well`: Integer well identifier (1-based index)
- `case`: Donor/case identifier (string)
- `name`: Read identifier from FASTQ
- `genomic_sequence`: Full read sequence
- `db_name`: Reference allele name (format: `GENE*ALLELE`)
- `gene`: Gene name (extracted from db_name by splitting on "*")
- `sequence`: Core allele sequence (may differ from db_seq for novel alleles)
- `count`: Number of supporting reads
- `ratio`: Allelic frequency (count / max_count_in_gene)
- `frequency`: Alternate name for ratio in some commands

## Search-Specific Columns

### exact
- `prefix`, `suffix`: Flanking sequences (length controlled by --affix)
- `heptamer`, `spacer`, `nonamer`: RSS elements (if --rss specified)
- `full_count`: Reads matching sequence + flanks
- `full_ratio`: Allelic ratio for full records
- `flank_index`: Variant number (1 to --top)
- `allele_freq`: Within-gene frequency
- `allele_case_freq`, `gene_case_freq`: Within-case frequencies
- `cross_case_median_*`: Cross-donor median values
- `*_to_cross_case_median_ratio`: Ratios for outlier detection
- `isin_db`: "" (known) or "Novel" (if --ref-fasta provided)

### blast
- `qseqid`, `sseqid`: BLAST query and subject IDs
- `qseq`: Original query sequence
- `db_seq`: Reference database sequence
- `aln_qseq`: Aligned and trimmed query sequence
- `aln_mismatch`: Mismatches in aligned portion
- `corecov`: Alignment coverage (length ratio)
- `isin_db`: Boolean, true if substring of known allele
- `allele_name`: Final assigned name (with _S hash if novel)
- BLAST columns: `pident`, `nident`, `length`, `mismatch`, `gapopen`, `qcovs`, `qcovhsp`, `qstart`, `qend`, `sstart`, `send`, `qlen`, `slen`, `evalue`, `bitscore`, `sstrand`

### hsmm
- `pre_nonamer`, `pre_spacer`, `pre_heptamer`: 5' RSS (9+12+7 bp)
- `post_heptamer`, `post_spacer`, `post_nonamer`: 3' RSS (7+12+9 bp)
- `heptamer_prob_pre`, `heptamer_prob_post`: PWM probabilities
- `heptamer_logp_pre`, `heptamer_logp_post`: Log-probabilities
- `log_path_prob`: Best path log-probability
- `log_total_prob`: Sum of all paths log-probability
- `posterior_prob`: Confidence (best path / total)
- `nearest_db`: Closest reference allele
- `nearest_db_dist`: Levenshtein distance
- `allele_name`: Final name (reference or _S hash)

## Analysis-Specific Columns

### association
- `allele_a`, `allele_b`: Allele pair
- `r`: Phi coefficient (Pearson correlation for binary presence/absence data across donors)
- `r2`: r squared (phi coefficient squared; equivalent to standardized LD r² but computed from unphased data)
- `jaccard`: Jaccard co-presence index
- `support`: n11 (co-present donors)
- `similarity`: Filtered similarity metric (r or r2 after applying jaccard/support filters)

### haplotype
- `genotype`: "homozygous", "heterozygous", "uncertain"
- `allele_1`, `allele_2`: Primary and secondary alleles
- `count_1`, `count_2`: Respective counts
- `ratio`: Minor/major allele ratio
- `total_count`: Total reads for gene in donor
- `other_alleles`: Additional alleles (if >2)
- `novel_1`, `novel_2`: Boolean flags (if --novel-fasta)

### bwa
- `position`: Chromosome and coordinate (format: "chr:pos")
- `edit_distance`: NM (mismatches in alignment)
- `ref_sequence`: Reference genome sequence (orientation-adjusted)
- `orientation`: "+" (forward) or "-" (reverse)

---

# Version History

This documentation reflects the current version of immunediscover. For version-specific changes, see the [releases page](https://gitlab.com/gkhlab/immunediscover.jl/-/releases).

**Deprecated Commands** (removed from current version):
- `hamming`: Hamming distance-based search
- `pattern`: Kmer-based pattern search
- `regex`: Regex-based insert detection
- `trim`: Profile-based trimming
- `nwpattern`: Needleman-Wunsch pattern alignment

Users of older versions should migrate to:
- `hamming` → `search exact` or `search blast`
- `pattern` → `search blast`
- `regex` → `search hsmm` (for D genes)
- `trim` → use trimming in `search exact` with --extension

---

# Support and Citation

For issues, bug reports, or feature requests, please use the [GitLab issue tracker](https://gitlab.com/gkhlab/immunediscover.jl/-/issues).

For questions about usage, please contact the development team or refer to the repository documentation.

If you use immunediscover in your research, please cite the repository and associated publications.
