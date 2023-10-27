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

## Standard
1. **Demultiplexing** of FASTQ file into correct cases
2. **Exact** Assign to demultiplexed reads exact matches to database of known alleles

## Advanced
Optional discovery steps can be performed to find novel alleles
3. **Pattern** search for alleles with kmers and trimming heuristics (faster but more false positives)
4. **Hamming** search similar sequences across all reads by sliding a window of known allele across the read (slower but more accurate)

## Demultiplexing
To demultiplex FASTQ file is is required to provide path to compressed FASTQ file (fastq.gz) and tab separated (TSV) file with following mandatory columns:
- forward_index
- reverse_index
- case
No other columns are accepted.

The following command separates reads to different cases by extracting forward_index on 5' end and reverse_index on 3' end using the same orientation as provided in `indices.tsv` file:

```
immunediscover demultiplex test.fastq.gz indices.tsv test.tsv
```

Output is a compressed tab separated (TSV) file with following columns:
- case
- name
- genomic_sequence

## Exact
This subcommand searches for exact matches to database of known alleles. It is required to provide path to compressed FASTQ file (fastq.gz) and tab separated (TSV) file with demultiplex data, the output is TSV file with the following columns:
- case
- db_name (name of the database sequence)
- count (number of reads matching the database sequence)
- gene
- frequency (count normalized by total number of reads in the case and for the gene)

To perform search for exact matches to database of known alleles:
```
immunediscover exact test.tsv.gz V.fasta test-exact.tsv.gz
```
## Hamming
This subcommand slides a window equal to database sequence length across entire read and performs assignment to closest match from the database based on **hamming distance**.
By default command operates on the `genomic_sequence` column which corresponds to full read, but process can be largerly accelerated by utilizing `trimmed_sequence` column which is added by `trim` command, this depends on the correct trimming though and might be less accurate, additionally with trimmed sequences no neigbouring heptamers can be extracted. 

Greedy hamming search
```
immunediscover hamming --maxdist 10 test.tsv.gz test.fasta test_hamming.tsv.gz
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

## Pattern
Can best be describe as an algorithm taking the following steps:
    - From starting database find kmers specific only to each input gene
    - Sample gene specific kmers
    - Search sequences containing these kmers
    - Group sequences by identity and apply filtering base on their count
    - Find levenshtein distance to known alleles and assign name based on the closest one

For example
```
immunediscover pattern test.tsv.gz V.fasta -l 200 -r 0.2 test-pattern.tsv.gz
```
Would perform analysis for entire plate, using trimmed genes of minimal length 200 and applying allelic ratio filter of 0.2.
Result is a compressed tab separated (TSV) file with candidate alleles and their counts for entire plate. Typically you want to pick some of these alleles and perform analysis on them separately by valindating them with hamming search or exact search module.

