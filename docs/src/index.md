```@meta
CurrentModule = immunediscover
```

# Immunediscover

[![Release](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mateusz-kaduk.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/commits/main)

## Description

Immunediscover is a package for immune repertoire sequencing. It can be downloaded as prebuild binary bundle from [releases](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases).

## Usage

Immunediscover has multiple commands that work as a pipeline and serve different purposes. The starting file the program is single FASTQ file with the reads and indexing barcodes (indices). Futher subcommands required demultiplex data in TSV format.

Typical application to find V genes takes following steps
1. **Demultiplexing** of FASTQ file into correct cases
2. **Heptamer** extraction and extension of provided query alleles

### Demultiplexing
To demultiplex FASTQ file is is required to provide path to compressed FASTQ file (fastq.gz) and tab separated (TSV) file with following mandatory columns:
- forward_index
- reverse_index
- case
No other columns are accepted.

The following command separates reads to different cases by extracting forward_index on 5' end and reverse_index on 3' end using the same orientation as provided in `indices.tsv` file:

```
immunediscover/bin/immunediscover demultiplex test.fastq.gz indices.tsv test.tsv
```

Output is a compressed tab separated (TSV) file with following columns:
- case
- name
- genomic_sequence

# API
Functions implemented in a package are documented and listed below:
```@index
```

```@autodocs
Modules = [immunediscover.demultiplex, immunediscover.data, immunediscover.simulate, immunediscover.profile, immunediscover.trim, immunediscover.cli]
```
