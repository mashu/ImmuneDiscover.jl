```@meta
CurrentModule = immunediscover
```

[![Release](https://gitlab.com/gkhlab/immunediscover.jl/-/badges/release.svg)](https://gitlab.com/gkhlab/immunediscover.jl/-/releases)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gkhlab.gitlab.io/immunediscover.jl/dev)
[![Build Status](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/pipeline.svg)](https://gitlab.com/gkhlab/immunediscover.jl/pipelines)
[![Coverage](https://gitlab.com/gkhlab/immunediscover.jl/badges/main/coverage.svg)](https://gitlab.com/gkhlab/immunediscover.jl/commits/main)

Immunediscover is a software package specifically developed for analyzing genomic next-generation sequencing (NGS) data of immune receptors. The package can be acquired as a prebuilt binary bundle from the [releases page](https://gitlab.com/gkhlab/immunediscover.jl/-/releases).

Immunediscover comprises multiple commands functioning as distinct standalone tools, each serving a unique purpose. The initial input for the program is a single FASTQ file containing the reads and TSV file with indexing barcodes (indices). Subsequent subcommands necessitate the demultiplexing of data into a TSV format.
# Installation
Downloading the Software: Immunediscover, a software package, is provided as a pre-compiled (prebuilt) binary. This means that you don't have to compile or build the software from its source code. You can download it directly from the releases page on their GitLab repository. Here is the link to that page:  [releases page](https://gitlab.com/gkhlab/immunediscover.jl/-/releases). Choose the version that matches your operating system.

1. Extracting the Downloaded File: After downloading, you will have a file that needs to be extracted (unzipped) to a folder on your computer. You can choose any directory you prefer. For example, if you choose to extract it to `/home/user/immunediscover`, you will find the program files in this directory.
2. Setting Up the Environment Variable:
    - The goal here is to make the immunediscover command accessible from any directory in your terminal. This is achieved by adding the path to the program's `bin` directory to an environment variable called `PATH`.
    - If you extracted the files to `/home/user/immunediscover`, the path to add will be `/home/user/immunediscover/bin/`.
    - You can temporarily update your `PATH` by running the following command in the terminal:
    ```
    export PATH=$PATH:/home/user/immunediscover/bin/
    ```
    - This command tells your system to add the `immunediscover/bin` directory to the list of places it checks for executable files.
3. Making the Change Permanent:
    - The export command you ran will only last for your current terminal session. If you open a new terminal window, it won't remember this change.
    - To make it permanent, you need to add this command to your `~/.bashrc` file, which is a script that runs every time you open a new terminal session.
    - You can edit this file with a text editor and add the export command to it. After saving and closing the file, the change will take effect every time you start a new terminal session.

# Usage
Basic usage typicallly consists of the following steps:
1. **Demultiplexing**: Segregates the FASTQ file into appropriate cases.
2. **Exact**: Assigns demultiplexed reads to exact matches within a database of known alleles.

These optional steps are aimed at discovering novel alleles:
3. **Pattern**: Conducts a search for alleles using kmers and trimming heuristics (this method is faster but may result in more false positives).
4. **Hamming**: Searches for similar sequences across all reads by sliding a window of a known allele across the read (this method does not discover differences in length but produces fewer false positives).

Some of the subcommands use multithreading to speed up the process. The number of threads can be controlled using the `JULIA_NUM_THREADS` environment variable. For instance, to use 4 threads, run the following command before executing immunediscover:
```bash
export JULIA_NUM_THREADS=4
```

# Demultiplexing
To demultiplex a FASTQ file, provide the path to the compressed FASTQ file (`fastq.gz`) and a tab-separated (TSV) file with the following mandatory columns:
- `forward_index`
- `reverse_index`
- `case`
No other columns are accepted.

The command below segregates reads into different cases based on the forward_index at the 5' end and reverse_index at the 3' end, aligning with the orientation provided in the `indices.tsv` file:

```bash
immunediscover demultiplex test.fastq.gz indices.tsv test.tsv
```

To run this command on multiple FASTQ files, consider naming libraries the index file the same way as the FASTQ file, e.g. `test.fastq.gz` and `test.tsv` for the first library, `test2.fastq.gz` and `test2.tsv` for the second library, etc. Then, use the following command to run the demultiplexing on all libraries:
```bash
for i in *.fastq.gz; do immunediscover demultiplex $i ${i%.fastq.gz}.tsv ${i%.fastq.gz}-demultiplex.tsv; done
```
Note that fastq files are assumed to be gzipped (which is good practice to save space). If they are not, use the following command:
```bash
for i in *.fastq; do immunediscover demultiplex $i ${i%.fastq}.tsv ${i%.fastq}-demultiplex.tsv; done
```
or altenatively, using findutils:
```bash
find . -name "*.fastq.gz" -exec bash -c 'immunediscover demultiplex $0 ${0%.fastq.gz}.tsv ${0%.fastq.gz}-demultiplex.tsv' {} \;
```

### Demultiplex Program Parameters

The `demultiplex` program uses the following parameters:

### Positional Arguments

1. `fastq`: FASTQ file with reads.
2. `indices`: TSV file with barcodes and cases.
3. `output`: TSV file to save demultiplex data.

### Optional Arguments

- `-l, --length LENGTH`: Minimum length of the read. (Default: 200)

- `-s, --split`: Split FASTQ files per case. This option does not have an associated type or default as it is a toggle switch. (Default: false)

### Demultiplex output
The output is a compressed TSV file with the following columns:
- `well`: The well on the plate.
- `case`: The case (donor) containing the allele.
- `name`: The name of the read.
- `genomic_sequence`: The sequence of the read.

# Exact
This subcommand identifies exact matches to a database of known alleles. Provide the path to the compressed FASTQ file (`fastq.gz`) and the TSV file with demultiplexed data. 

To execute a search for exact matches to a database of known alleles, use the command:
```bash
immunediscover exact -g V -c 10 -f 0.1 test.tsv.gz V.fasta test-exact.tsv.gz
```
Above command searches for exact matches to the database of V alleles, requiring a minimum of `10` reads per allele and a minimum frequency of `10%`. The result is a compressed TSV file containing exact matches for alleles and their counts for each well and case.

### Exact Program Parameters

The `exact` program uses the following parameters:

### Positional Arguments

1. `tsv`: TSV file with demultiplexed data.
2. `fasta`: FASTA file with query alleles.
3. `output`: TSV file to save output.

The records stored in the output can be categorized as follows:
- Collapsed Records: These are matches to the database that exclude the flanking regions.
- Full Records: These include matches to the database along with the flanking regions.

Generally, the counts for full records tend to be smaller from those of collapsed records. This difference arises because flanking regions are often more unique and dilute the total counts as compared to records collapsed by the allelle sequence only.

### Optional Arguments
- `-f, --minfreq FREQ`: Minimum allelic ratio applied within each gene group. (Default: 0.01)
- `-c, --mincount MINCOUNT`: Minimum cluster size. (Default: 10)
- `--full-minfreq FREQ`: Minimum allelic ratio applied within each gene group for **full record**. (Default: 0.01)
- `--full-mincount MINCOUNT`: Minimum cluster size for **full record**. (Default: 1)
- `-n, --noplot`: Disable unicode gene plot. This option does not have associated types or defaults as it is a toggle switch.
- `-g, --gene`: Gene to use for plotting. (Default: "V")
- `-u, --uncollapsed`: Save uncollapsed full records. (Default: false)

### Exact output
The output is a TSV file containing the following columns:
- `well`: Identifies the well on the plate.
- `case`: Specifies the case (donor) that contains the allele.
- `db_name`: Denotes the name of the allele as listed in the database.
- `sequence`: Represents the sequence of the allele.
- `full_count`: Indicates the number of full records in which the allele matched the read.
- `count`: Reflects the number of records in which the allele matched the read.
- `gene`: Specifies the gene associated with the allele.
- `full_frequency`: Describes the frequency of full records in which the allele matched the read.
- `frequency`: Captures the frequency of records in which the allele matched the read.
- `log2_count`: Represents the logarithm base 2 of the count.

Additional columns, such as `prefix`, `suffix`, `heptamer`, `spacer`, and `nonamer`, are included based on the `--gene` option, which is set to `V` by default.

# Hamming

This subcommand executes a sliding window operation, equal to the database sequence length, across the entire read to find the closest match from the database based on hamming distance. By default, the command operates on the `genomic_sequence` column corresponding to the full read. To expedite the process, utilize the `trimmed_sequence` column added by the trim command, although this may compromise accuracy and eliminate the extraction of neighboring heptamers.

For a greedy Hamming search, use the command:
```bash
immunediscover hamming --maxdist 10 test.tsv.gz test.fasta test_hamming.tsv.gz
```

### Hamming Program Parameters

The `hamming` program accepts the following parameters:

### Positional Arguments

1. `tsv`: TSV file with demultiplexed data.
2. `fasta`: FASTA file with query alleles.
3. `output`: TSV file to save output.

### Optional Arguments

- `-a, --assignments ASSIGNMENTS`: Optional TSV file to save intermediate assignments.
- `-d, --maxdist MAXDIST`: Maximum distance allowed. (Default: 2)
- `-c, --mincount MINCOUNT`: Minimum cluster size. (Default: 10)
- `-r, --ratio RATIO`: Minimum allelic ratio applied on cluster sizes. (Default: 0.25)
- `-f, --column COLUMN`: Column with genomic sequence. (Default: "genomic_sequence")
- `-u, --umi`: Indicates that UMI (Unique Molecular Identifier) is present in the read. (Default: false)
- `-l, --limit LIMIT`: Limit to this number of sequences, zero means no limit. (Default: 0)
- `-n, --noplot`: Disable unicode gene plot. This option does not have an associated type or default as it is a toggle switch. (Default: false)

### Hamming output
The output is a compressed TSV file with the following columns:
- `gene`: The gene of the allele.
- `prefix`: The prefix sequence of the allele.
- `middle`: The sequence of the allele.
- `suffix`: The suffix sequence of the allele.
- `count`: The number of reads matching the allele.
- `allele_name`: The new name of the allele.
- `case`: The case (donor) containing the allele.
- `well`: The well on the plate.

Note that `count` refers to the combined occurrences of `prefix`, `middle`, and `suffix` having the middle part best matching the database sequence.

# Pattern
![Pattern-Diagram](assets/pattern-diagram.png)

This subcommand can be elucidated through the following algorithmic steps:
1. Identify gene-specific kmers from the starting database that are unique to each gene but not part of any pseudo-gene.
2. Sample gene-specific kmers.
3. Search for sequences containing sampled kmers.
4. Trim matching reads to the gene sequence using the profile at the 5' and 3' ends.
5. Group sequences by identity and apply count-based filtering.
6. Calculate the Levenshtein distance to known alleles and assign names based on the closest match.

For instance, the command:
```bash
immunediscover pattern test.tsv.gz V.fasta -b pseudogenes.fasta -l 200 -f 0.1 -c 25 -w 30 -d 10 test-pattern.tsv.gz
```
will analyze the entire plate using with the following specifications:
- `-l 200`: Sets the minimum gene sequence length for consideration to 200 nucleotides.
- `-f 0.1`: Filters alleles, keeping those with a frequency of at least 10%.
- `-c 25`: Includes only alleles observed a minimum of 25 times.
- `-w 30`: Length of the profile on both 5' and 3' ends of a gene used for trimming reads to the gene sequence, set at 30.
- `-d 10`: Maximum distance for an allele to be considered novel, set at 10. Helps to exclude vastly different sequences like pseudogenes and spurious matches.
The result is a compressed TSV file containing candidate alleles and their counts for the entire plate. Typically, select some of these alleles for further analysis using Exact search.

### Pattern Program Parameters

The `pattern` program accepts the following parameters:

### Positional Arguments

1. `input`: TSV file with demultiplex data.
2. `fasta`: FASTA file with aligned gene sequences to build trimming profile.
3. `output`: TSV file to save demultiplex data.

### Optional Arguments

- `-t, --top TOP`: Save top candidates with highest counts per gene. (Default: 5)
- `-w, --weights WEIGHTS`: Length of the position weight matrix. (Default: 20)
- `-n, --noprofile`: Use distance to database gene lengths instead of profiles to trim reads. (Default: false)
- `-l, --length LENGTH`: Minimum length of the trimmed read. (Default: 200)
- `-k, --kmer KMER`: Kmer size which will be used to search for patterns. (Default: 12)
- `-m, --maxkmer MAXKMER`: Maximum kmer size if kmer size needs to be increased automatically. (Default: 50)
- `-d, --maxdist MAXDIST`: Maximum distance allowed for alleles. (Default: 50)
- `-s, --sample SAMPLE`: Number of kmers to sample from the set of all kmers to search for a gene. (Default: 5)
- `-f, --minfreq MINFREQ`: Minimum allelic ratio applied within each gene group. (Default: 0.01)
- `-c, --mincount MINCOUNT`: Minimum count for an allele. (Default: 10)
- `-b, --blacklist BLACKLIST`: Blacklist file with sequences to be excluded from pattern search (e.g. pseudo-genes).

### Pattern output
The output is a compressed TSV file with the following columns:
- `count`: The number of reads matching the allele.
- `closest`: The closest allele from the database.
- `dist`: The Levenshtein distance to the closest allele.
- `allele_name`: The new name of the allele.
- `gene`: The gene of the allele.
- `length`: The length of the allele.
- `length_db`: The length of the closest allele from the database.
- `patterns`: The patterns used to find the allele.
- `ratio`: The ratio of the allele with respect to highest allele in the gene group.
- `seq`: The sequence of the allele.

# Regex
The Regex search command operates on a FASTA format database to find exact matches in DNA sequences, extracting the prefixes and suffixes of matching alleles. It then generates an array of all prefix and suffix combinations, searching within the sequences to identify inserts between flanking sequences. The result is a TSV file compressed to include potential alleles, along with their occurrence counts for each well and case. This method is particularly effective in identifying short, unknown alleles such as D genes. Subsequently, the identified insert sequences are aligned using the Smith-Waterman algorithm against a known alleles database. This step assigns names to the closest matches and performs both consolidation and filtering. The final output provides a comprehensive list of alleles corresponding to the genotype of each well and case.

To analyze the plate using regex search, use the command:
```bash
immunediscover regex test.tsv.gz test.fasta test_regex.tsv.gz
```

### Positional Arguments

1. `tsv`: TSV file with demultiplexed data.
2. `fasta`: FASTA file with query alleles.
3. `output`: TSV file to save output.

### Optional Arguments
- `--insert-minlen`: Minimum length of an insert (Default: 15)
- `--insert-maxlen`: Maximum length of an insert (Default: 40)
- `--flank-mincount`: Minimum number of reads per flanks (Default: 25)
- `--flank-frequency`: Minimum frequency of flanks (Default: 0.5)
- `-c, --mincount`: Minimum count of a match (Default: 3)
- `-f,--frequency`: Lowest frequency of a match (Default: 0.1)
- `-p, --nprefix`: Number of nucleotides to extract from the 5' end of the query sequence (Default: 7)
- `-s, --nsuffix`: Number of nucleotides to extract from the 3' end of the query sequence (Default: 7)

### Regex output
The output is a compressed TSV file with the following columns:
- `well`: The well on the plate.
- `case`: The case (donor) containing the allele.
- `best_name`: The name of the allele.
- `distance`: The distance to the closest allele from the database.
- `prefix`: The prefix sequence of the allele.
- `best_aln`: The alignment of the allele to the database.
- `suffix`: The suffix sequence of the allele.
- `count`: The number of reads matching the allele.
- `gene`: The gene of the allele.
- `frequency`: The frequency of the allele.

# Collect
The Collect command is used to collect the results of multiple runs of the Exact, Hamming, Pattern, and Regex commands. The command accepts a list of TSV files and outputs a single TSV file with the combined results.

To collect the results of multiple runs, use the command:
```bash
immunediscover collect test-exact*.tsv.gz exact-collected.tsv
```

### Positional Arguments

1. `pattern`: Pattern of TSV files to collect (i.e. `test-exact*.tsv.gz`, `test-regex*.tsv.gz`, `test-pattern*.tsv.gz`...).
2. `output`: Output TSV file to save collected data.

### Collect output
The output is a compressed TSV file with the same columns as the input files and an additional column `file` indicating the source of the record.

# Bwa
The Bwa command is used to align reads to a reference genome using the BWA-MEM algorithm. The command accepts a TSV file with `allele_name` and `sequence` columns and outputs filtered TSV file for sequences that match requested `chromosome`.

To align reads to a reference genome, use the command:
```bash
bwa index gnome.fasta  # This step is required only once per reference genome
immunediscover bwa test-regex.tsv.gz gnome.fasta test_bwa.tsv.gz -c "chromosome 14"
```

### Positional Arguments
1. `tsv`: TSV file with columns allele_name and seq
2. `genome`: FASTA file with indexed genome
3. `output`: TSV file to save filtered input

### Optional Arguments
- `-c, --chromosome`: Chromosome to align to that has this string in genomic record description (Default: "chromosome 14")
- `-n, --colname`: Column name to use for allele name (Default: "allele_name")
- `-s, --colseq`: Column names to use for sequence (Default: "prefix best_aln suffix")
