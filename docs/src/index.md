```@meta
CurrentModule = immunediscover
```

# Immunediscover

[![Release](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mateusz-kaduk.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/commits/main)

Immunediscover is a package for immune repertoire sequencing. It can be downloaded as prebuild binary bundle from [releases](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases).

# Usage

Immunediscover has multiple commands that work as separate stand-alone tools and serve different purposes. The starting file the program is single FASTQ file with the reads and indexing barcodes (indices). Futher subcommands required demultiplex data in TSV format.

Typical application to find V genes takes following steps
1. **Demultiplexing** of FASTQ file into correct cases
2. **Pattern** search for alleles that may contain inserts or indels and call novel alleles
3. **Hamming** search similar sequences across all reads (not allowing for indels or inserts) and call novel alleles
4. **Exact** Fast tool to search demultiplexed data for exact matches to database of known alleles collected from previous steps

## Demultiplexing
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

## Trim (optional)
To trim demultiplex data in TSV format with `genomic_sequence` column are needed as input as well as FASTA file with V alleles.
Program iterates reads and attempts to find start and stop position of potential V allele in each read, reports fraction of succesfully trimmed sequences (for start < stop) and appends additional column `trimmed_sequence`. To save space, it is possible to save `start` and `stop` positions instead with a subcommand switch. 

```
immunediscover/bin/immunediscover trim test.tsv.gz test.fasta test_trimmed.tsv.gz -l 6
```

**Note** for V genes longer profile `-l 20` worked well with our libraries, example above uses shorter profile for generated test sequences.

Output is a compressed tab separated (TSV) file with following columns:
- case
- name
- genomic_sequence
- trimmed_sequence

## Pattern
Can best be describe as an algorithm taking the following steps:
    - From starting database find kmers specific only to each input gene
    - Sample gene specific kmers
    - Search sequences containing these kmers
    - Group sequences by identity and apply filtering base on their count
    - Find levenshtein distance to known alleles and assign name based on the closest one

For example
```
immunediscover pattern COL.tsv.gz V_Immune24Oct2023.fasta -l 200 -r 0.2 COL-pattern.tsv.gz
```
Would perform analysis for entire plate, using trimmed genes of minimal length 200 and applying allelic ratio filter of 0.2.
Result is a compressed tab separated (TSV) file with candidate alleles and their counts for entire plate. Typically you want to pick some of these alleles and perform analysis on them separately by valindating them with hamming search or exact search module.

## Hamming
This subcommand slides a window equal to database sequence length across entire read and performs assignment to closest match from the database based on **hamming distance**.
By default command operates on the `genomic_sequence` column which corresponds to full read, but process can be largerly accelerated by utilizing `trimmed_sequence` column from the previous step. 

Faster search using trimming heuristic
```
immunediscover/bin/immunediscover hamming --maxdist 10 -f trimmed_sequence test_trimmed.tsv.gz test.fasta test_hamming.tsv.gz
```

Slower but greedy search which is more accurate
```
immunediscover/bin/immunediscover hamming --maxdist 10 test_trimmed.tsv.gz test.fasta test_hamming.tsv.gz
```

Output is a compressed tab separated (TSV) file with following columns:
- gene
- prefix
- middle
- suffix
- cluster_size
- allele_name
- case
Where `cluster_size` is the number of combined `prefix`,`middle` and `suffix` occurences that have `middle` part best matching database sequence.
