module cli
    using ArgParse
    export parse_commandline
    
    """
        always_gz(file_path)

    Return path that ends with .gz or with appeneded extension
    """
    function always_gz(file_path)
        if endswith(file_path, ".gz")
            output = file_path
        else
            output = file_path*".gz"
        end
    end

    """
        parse_commandline(args)

    Handle command line
    """
    function parse_commandline(args)
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "0.1",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "GKHLab, 2023")
        @add_arg_table! s begin
            "demultiplex"
                help = "Demultiplex plate library"
                action = :command
            "heptamer"
                help = "Trim and extend Vs up to heptamer"
                action = :command
            "trim"
                help = "Trim reads with start and stop profiles"
                action = :command
            "exact"
                help = "Search for exact matches"
                action = :command
            "hamming"
                help = "Search for matches within hamming distance"
                action = :command
            "pattern"
                help = "Search for novel alleles with kmers and trim with PWM"
                action = :command
            end

        @add_arg_table! s["demultiplex"] begin
            "fastq"
                help = "FASTQ file with reads"
                required = true
            "indices"
                help = "TSV file with barcodes and cases"
                required = true
            "output"
                help = "TSV file to save demultiplex data"
                required = true
            "-l", "--length"
                help = "Minimum length of the read"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 200
            "-s", "--split"
                help = "Split FASTQ files per case"
                action = :store_true
        end

        choices = ["IGKV", "IGLV", "IGHV"]
        @add_arg_table! s["heptamer"] begin
            "tsv"
                help = "TSV file with demultiplexed reads"
                required = true
            "fasta"
                help = "FASTA file with query alleles"
                required = true
            "output"
                help = "TSV file to save data with identified heptamers"
                required = true
            "summary"
                help = "TSV file to save summary collapsed alleles with statistics"
                required = true
            "-j", "--json"
                help = "JSON file with dictionary containing haptamers"
                default = "heptamers.json"
            "-c", "--chain"
                default = "IGKV"
                range_tester = (x->x ∈ choices)
                arg_type = String
                help = "chain; must be one of " * join(choices, ", ", " or ")
            "-d", "--maxdist"
                help = "A positive integer indicating maximum Hamming distance from any of heptamers in JSON file"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 1
            "-b", "--begin"
                help = "How much to trim from the 5' beginning of the query sequence"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 0
            "-e", "--end"
                help = "How much to trim from the 3' end of the query sequence"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 8
            "-m", "--mincount"
                help = "Minimum count allowed in summary"
                default = 1
                arg_type = Int
                range_tester = (x->x >= 1)
            "-r","--ratio"
                help = "Lowest allowed ratio between counts of full allele sequence and trimmed allele sequence"
                default = 0.25
                arg_type = Float64
                range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        @add_arg_table! s["hamming"] begin
        "tsv"
            help = "TSV file with demultiplexed data"
            required = true
        "fasta"
            help = "FASTA file with query alleles"
            required = true
        "output"
            help = "TSV file to save ouput"
            required = true
        "-a", "--assignments"
            help = "Optional TSV file to save intermediate assignments"
            arg_type = String
        "-d", "--maxdist"
            help = "Maximum distance allowed"
            default = 2
            arg_type = Int
            range_tester = (x->x >= 0)
        "-c", "--mincount"
            help = "Minimum cluster size"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "-r","--ratio"
            help = "Minimum allelic ratio applied on cluster sizes"
            default = 0.25
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-f","--column"
            help = "Column with genomic sequence"
            default = "genomic_sequence"
            arg_type = String
        "-u", "--umi"
            help = "UMI is present in the read"
            action = :store_true
        "-l", "--limit"
            help = "Limit to this number of sequences, zero means no limit"
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "-n","--noplot"
            help = "Disable unicode gene plot"
            action = :store_true
        end

        @add_arg_table! s["trim"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with aligned gene sequences to build trimming profile"
            required = true
        "output"
            help = "TSV file to save demultiplex data"
            required = true
        "-w", "--weights"
            help = "Length of the position weight matrix"
            default = 20
            arg_type = Int
            range_tester = (x->x > 1)
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 1
            arg_type = Int
            range_tester = (x->x >= 1)
        "-p", "--position"
            help = "An optional switch to save start,stop columns instead of trimmed_sequence"
            action = :store_true
        end

        @add_arg_table! s["exact"] begin
        "tsv"
            help = "TSV file with demultiplexed data"
            required = true
        "fasta"
            help = "FASTA file with query alleles"
            required = true
        "output"
            help = "TSV file to save ouput"
            required = true
        "-r","--ratio"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.01
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c", "--mincount"
            help = "Minimum cluster size"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "-n","--noplot"
            help = "Disable unicode gene plot"
            action = :store_true
        end

        @add_arg_table! s["pattern"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with aligned gene sequences to build trimming profile"
            required = true
        "output"
            help = "TSV file to save demultiplex data"
            required = true
        "-w", "--weights"
            help = "Length of the position weight matrix"
            default = 20
            arg_type = Int
            range_tester = (x->x > 1)
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 1
            arg_type = Int
            range_tester = (x->x >= 1)
        "-k", "--kmer"
            help = "Kmer size which will be used to search for patterns"
            default = 12
            arg_type = Int
            range_tester = (x->x >= 3)
        "-m", "--maxkmer"
            help = "Maximum kmer size if kmer size needs to be increased automatically"
            default = 50
            arg_type = Int
            range_tester = (x->x >= 3)
        "-s", "--sample"
            help = "Number of kmers to sample from the set of all kmers to search for a gene"
            default = 4
            arg_type = Int
            range_tester = (x->x >= 1)
        "-f","--minfreq"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.01
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c","--mincount"
            help = "Minimum count for an allele"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "-b","--blacklist"
            help = "Blacklist file with sequences to be excluded from pattern search (e.g. pseudo-genes)"
            arg_type = String
        end
        # Show help if no arguments are provided
        try
            return parse_args(args, s)
        catch e
            if isa(e, ArgParseError)
                println(e)
                ArgParse.show_help(s)
            else
                rethrow()
            end
        end
    end
end
