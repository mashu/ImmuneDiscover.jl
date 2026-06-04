# Immunediscover

[![Release](https://gitlab.com/gkhlab/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/gkhlab/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gkhlab.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/gkhlab/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/gkhlab/immunediscover.jl/commits/main)
[![Codecov](https://codecov.io/gl/gkhlab/immunediscover.jl/graph/badge.svg?token=K5GWGJP8FS)](https://codecov.io/gl/gkhlab/immunediscover.jl)

Immunediscover analyzes genomic NGS data from immune receptors to discover and
characterize immunoglobulin alleles (V/D/J). Download a prebuilt binary from the
[releases](https://gitlab.com/gkhlab/immunediscover.jl/-/releases) page, or build
from source (see below).

## Command groups

| Group | Purpose | Subcommands |
|-------|---------|-------------|
| `preprocess` | Barcode demultiplexing | `demultiplex` |
| `discover` | De novo allele discovery | `blast`, `hsmm` |
| `search` | Search against known references | `exact`, `heptamer`, `bwa` |
| `analyze` | Downstream analysis | `cooccurrence`, `haplotype` |
| `table` | TSV utilities | `outerjoin`, `leftjoin`, `transform`, `aggregate`, `unique`, `sort`, `filter`, `select`, `fasta`, `collect`, `exclude` |
| `fasta` | FASTA utilities | `merge`, `diff`, `hash` |

## Quick start

```bash
# 1. Demultiplex reads by plate barcodes
immunediscover preprocess demultiplex reads.fastq.gz indices.tsv demux.tsv.gz

# 2. Exact-match reads to a known reference (V genes)
immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V

# 3. Discover novel V alleles with BLAST (gene preset)
immunediscover discover blast demux.tsv.gz IGHV.fasta blast_V.tsv.gz -g V

# 4. Detect (short) D genes with the RSS HSMM
immunediscover discover hsmm demux.tsv.gz IGHD.fasta hsmm_D.tsv.gz

# 5. Export discovered alleles to FASTA
immunediscover table fasta blast_V.tsv.gz novel_V.fasta --filter Novel
```

`immunediscover <group> <subcommand> --help` documents every option. Full
documentation lives under [`docs/`](docs/src) and the
[dev docs site](https://gkhlab.gitlab.io/immunediscover.jl/dev).

## From source

```bash
# Run without building (uses the project environment)
./scripts/run.sh preprocess demultiplex reads.fastq.gz indices.tsv demux.tsv.gz

# Run the test suite
julia --project -e 'using Pkg; Pkg.test()'

# Build a standalone binary (PackageCompiler)
julia --project=build scripts/build_binary.jl
```
