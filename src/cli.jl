module cli
    using ArgParse

    """
    Return path that ends with .gz or with appeneded extension
    """
    function always_gz(file_path)
        if endswith(file_path, ".gz")
            output = file_path
        else
            output = file_path*".gz"
        end
    end

    export parse_commandline
    """
    Handle command line
    """
    function parse_commandline(args)
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "0.1",
                            add_version = true,
                            epilog = "Mateusz Kaduk <mateusz.kaduk@ki.se>")
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
        "-l", "--length"
            help = "Profile length"
            default = 20
            arg_type = Int
            range_tester = (x->x > 1)
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
        end

        return parse_args(args,s)
    end
end