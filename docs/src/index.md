```@meta
CurrentModule = immunediscover
```

[![Release](https://gitlab.com/gkhlab/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/gkhlab/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gkhlab.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/gkhlab/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/gkhlab/immunediscover.jl/commits/main)

# Immunediscover

Immunediscover is a comprehensive software package for analyzing genomic next-generation sequencing (NGS) data from immune receptors, specifically designed for discovering and characterizing immunoglobulin alleles.

## Overview

### Command Groups

Immunediscover provides a complete pipeline organized into command groups:

| Group | Purpose | Subcommands |
|-------|---------|-------------|
| **preprocess** | Barcode-based read segregation | demultiplex |
| **discover** | De novo allele discovery | blast, hsmm |
| **search** | Search against known references | exact, heptamer, bwa |
| **analyze** | Downstream analysis | cooccurrence, haplotype |
| **table** | TSV data manipulation | outerjoin, leftjoin, transform, aggregate, unique, sort, filter, select, fasta, collect, exclude |
| **fasta** | FASTA file operations | merge, diff, hash |

### Typical Pipeline

```
FASTQ + indices → preprocess → discover/search → analyze → export
```

1. **preprocess demultiplex** raw reads by plate barcodes
2. **discover** novel alleles (blast, hsmm) or **search** known references (exact, heptamer)
3. **analyze** results (co-occurrence, haplotypes) and QC with **search bwa**
4. **Export** to FASTA (`table fasta` / `fasta merge`) for validation

---

## Installation

### Download and Setup

1. **Download** pre-compiled binary from [releases page](https://gitlab.com/gkhlab/immunediscover.jl/-/releases)

2. **Extract** to directory (e.g., `/home/user/immunediscover`)

3. **Add to PATH**:
   ```bash
   export PATH=$PATH:/home/user/immunediscover/bin/
   echo 'export PATH=$PATH:/home/user/immunediscover/bin/' >> ~/.bashrc
   ```

4. **Verify**:
   ```bash
   immunediscover --version
   ```

### Performance Settings

Enable multithreading for faster processing:
```bash
export JULIA_NUM_THREADS=16  # Use available CPU cores
```

---

## Quick Start

### Minimal Example

```bash
# 1. Demultiplex by barcodes
immunediscover preprocess demultiplex reads.fastq.gz indices.tsv demux.tsv.gz

# 2. Exact match to known V genes
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V

# 3. Discover novel V genes
immunediscover discover blast demux.tsv.gz IGHV.fasta blast_V.tsv.gz -g V

# 4. Export novel alleles
immunediscover table fasta blast_V.tsv.gz novel_V.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 2
```

### Required Inputs

- **FASTQ**: Single-end reads (`.fastq.gz` or `.fastq`)
- **Indices TSV**: Barcodes with columns `forward_index`, `reverse_index`, `case`
- **Database FASTA**: Reference alleles (e.g., IMGT database)

---

## Next Steps

- See [Commands](commands.md) for detailed command reference
- See [Workflows](workflows.md) for complete analysis pipelines
- See [Parameter Guidelines](parameters.md) for threshold selection
- See [Troubleshooting](troubleshooting.md) for common issues
