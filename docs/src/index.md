```@meta
CurrentModule = immunediscover
```

# Immunediscover
[![Release](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mateusz-kaduk.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/mateusz-kaduk/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/mateusz-kaduk/immunediscover.jl/commits/main)

ImmuneDiscover is a software package designed for immune repertoire sequencing analysis. The package can be acquired as a prebuilt binary bundle from the [releases page](https://gitlab.com/mateusz-kaduk/immunediscover.jl/-/releases).

# Usage
ImmuneDiscover comprises multiple commands functioning as distinct standalone tools, each serving a unique purpose. The initial input for the program is a single FASTQ file containing the reads and indexing barcodes (indices). Subsequent subcommands necessitate the demultiplexing of data into a TSV format.

## Standard Operations
1. **Demultiplexing**: Segregates the FASTQ file into appropriate cases.
2. **Exact**: Assigns demultiplexed reads to exact matches within a database of known alleles.

## Advanced Operations
These optional steps are aimed at discovering novel alleles:
3. **Pattern**: Conducts a search for alleles using kmers and trimming heuristics (this method is faster but may result in more false positives).
4. **Hamming**: Searches for similar sequences across all reads by sliding a window of a known allele across the read (this method is slower but more precise).

## Demultiplexing
To demultiplex a FASTQ file, provide the path to the compressed FASTQ file (`fastq.gz`) and a tab-separated (TSV) file with the following mandatory columns:
- forward_index
- reverse_index
- case
No other columns are accepted.

The command below segregates reads into different cases based on the forward_index at the 5' end and reverse_index at the 3' end, aligning with the orientation provided in the `indices.tsv` file:

```bash
immunediscover demultiplex test.fastq.gz indices.tsv test.tsv
```
The output is a compressed TSV file with the following columns:
- case
- name
- genomic_sequence

## Exact
This subcommand identifies exact matches to a database of known alleles. Provide the path to the compressed FASTQ file (`fastq.gz`) and the TSV file with demultiplexed data. The output is a TSV file with the following columns:
- case
- db_name (name of the database sequence)
- count (number of reads matching the database sequence)
- gene
- frequency (count normalized by the total number of reads in the case and for the gene)

To execute a search for exact matches to a database of known alleles, use the command:
```bash
immunediscover exact test.tsv.gz V.fasta test-exact.tsv.gz
```

## Hamming

This subcommand executes a sliding window operation, equal to the database sequence length, across the entire read to find the closest match from the database based on hamming distance. By default, the command operates on the `genomic_sequence` column corresponding to the full read. To expedite the process, utilize the `trimmed_sequence` column added by the trim command, although this may compromise accuracy and eliminate the extraction of neighboring heptamers.

For a greedy Hamming search, use the command:
```bash
immunediscover hamming --maxdist 10 test.tsv.gz test.fasta test_hamming.tsv.gz
```
The output is a compressed TSV file with the following columns:
- gene
- prefix
- middle
- suffix
- cluster_size
- allele_name
- case

Here, `cluster_size` refers to the combined occurrences of `prefix`, `middle`, and `suffix` having the middle part best matching the database sequence.

## Pattern
This subcommand can be elucidated through the following algorithmic steps:
- Identify gene-specific kmers from the starting database.
- Sample these gene-specific kmers.
- Search for sequences containing these kmers.
- Group sequences by identity and apply count-based filtering.
- Calculate the Levenshtein distance to known alleles and assign names based on the closest match.

For instance, the command:
```bash
immunediscover pattern test.tsv.gz V.fasta -b pseudogenes.fasta -l 200 -r 0.2 test-pattern.tsv.gz
```
will analyze the entire plate using trimmed genes of a minimum length of 200, applying an allelic ratio filter of 0.2. The result is a compressed TSV file containing candidate alleles and their counts for the entire plate. Typically, select some of these alleles for further analysis using the Hamming or Exact search module.
