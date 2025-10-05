module cli
    using ArgParse
    export parse_commandline, apply_blast_presets!, show_blast_presets, show_blast_params
    export BLAST_PRESETS, BLAST_CLI_DEFAULTS
    import ArgParse.parse_item
    using ArgParse: @add_arg_table!
    using Logging
    using Dates

    # CLI defaults for blast command - single source of truth (based on v0.0.66)
    const BLAST_CLI_DEFAULTS = Dict(
        "forward" => 20,
        "reverse" => 20,
        "minfullfreq" => 0.1,
        "length" => 290,
        "maxdist" => 20,
        "minfullcount" => 5,
        "edge" => 0,
        "subjectcov" => 0.1,
        "minquality" => 0.75,
        "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
    )

    const BLAST_PRESETS = Dict(
        "V" => Dict(
            "forward" => 12,
            "reverse" => 12,
            "minfullfreq" => 0.1,
            "length" => 290,
            "maxdist" => 10,
            "minfullcount" => 10,
            "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
        ),
        "D" => Dict(
            "forward" => 40,
            "reverse" => 40,
            "minfullfreq" => 0.2,
            "length" => 5,
            "maxdist" => 20,
            "minfullcount" => 10,
            "edge" => 10,
            "subjectcov" => 0.25,
            "minquality" => 0.5,
            "args" => "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5"
        ),
        # Old D presets (commented out for reference)
        # "D" => Dict(
        #     "forward" => 0,
        #     "reverse" => 0,
        #     "minfullfreq" => 0.2,
        #     "length" => 10,
        #     "maxdist" => 10,
        #     "minfullcount" => 10,
        #     "args" => "-task blastn-short -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"
        #     #"args" => "-task blastn -word_size 7 -evalue 100 -penalty -2 -reward 1 -dust no -soft_masking false -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"
        # ),
        "J" => Dict(
            "forward" => 12,
            "reverse" => 12,
            "minfullfreq" => 0.1,
            "length" => 10,
            "maxdist" => 10,
            "minfullcount" => 10,
            "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"
        )
    )

    function show_blast_presets()
        println("BLAST Presets:")
        for (gene, settings) in BLAST_PRESETS
            println("\n$gene gene settings:")
            for (param, value) in settings
                println("  --$param = $value")
            end
        end
    end

    function show_blast_params(args)
        if args["%COMMAND%"] == "blast"
            gene = args["blast"]["gene"]
            if haskey(BLAST_PRESETS, gene)
                preset = BLAST_PRESETS[gene]
                println("BLAST Preset for $gene gene:")
                for (param, value) in preset
                    println("  --$param = $(args["blast"][param])")
                end
            end
        end
    end



    function apply_blast_presets!(parsed_args)
        if parsed_args["%COMMAND%"] == "blast"
            gene = parsed_args["blast"]["gene"]
            if haskey(BLAST_PRESETS, gene)
                preset = BLAST_PRESETS[gene]

                # Only apply preset if the current value equals the CLI default
                # This means the user didn't explicitly override it
                for (key, preset_value) in preset
                    if haskey(parsed_args["blast"], key) && haskey(BLAST_CLI_DEFAULTS, key)
                        current_value = parsed_args["blast"][key]
                        default_value = BLAST_CLI_DEFAULTS[key]
                        
                        # Apply preset only if current value equals the default
                        if current_value == default_value
                            @info "Applying $gene preset $key: $current_value → $preset_value"
                            parsed_args["blast"][key] = preset_value
                        else
                            @info "Keeping user override for $key: $current_value (not default $default_value)"
                        end
                    else
                        # Apply preset for keys not in CLI defaults
                        if haskey(parsed_args["blast"], key)
                            @info "Applying $gene preset $key: $(parsed_args["blast"][key]) → $preset_value"
                            parsed_args["blast"][key] = preset_value
                        end
                    end
                end
            end
        end
        return parsed_args
    end

    function get_latest_git_tag()
	    return strip(read(`git -C $(@__DIR__) describe --tags --abbrev=0 `, String))
    end

    const SOFTWARE_VERSION = get_latest_git_tag()
    const SOFTWARE_GIT_HASH = strip(read(`git -C $(@__DIR__) rev-parse HEAD`, String))

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
                            version = "$SOFTWARE_VERSION (git $SOFTWARE_GIT_HASH)",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "GKHLab, $SOFTWARE_VERSION (git $SOFTWARE_GIT_HASH)")
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
            "nwpattern"
                help = "Search for novel alleles with nwpattern and trim with PWM"
                action = :command
            "pattern"
                help = "Search for novel alleles with kmers and trim with PWM"
                action = :command
            "regex"
                help = "Search for short novel alleles with regex"
                action = :command
            "hsmm"
                help = "Detect D genes using an HSMM trained on RSS flanks"
                action = :command
            "fasta"
                help = "Extract sequences from TSV file into FASTA format"
                action = :command
            "merge"
                help = "Merge multiple FASTA files, keeping only unique sequences"
                action = :command
            "hash"
                help = "Add hash based _S suffix to all allele names in the FASTA file"
                action = :command
            "collect"
                help = "Helper function to collect all TSV files into one"
                action = :command
            "exclude"
                help = "Exclude sequences from the TSV file"
                action = :command
            "bwa"
                help = "Filter input TSV file by mapping to the genome"
                action = :command
            "blast"
                help = "Run BLAST on the input TSV file and perform identity clustering for discovery"
                action = :command
            "diff"
                help = "Diff two FASTA files based on sequence identity but keep associated names"
                action = :command
            "linkage"
                help = "Analyze cross-donor co-occurrence (linkage) between alleles"
                action = :command
            "haplotype"
                help = "Infer approximate haplotypes from unphased data based on diploid assumptions"
                action = :command
            "table"
                help = "Table operations"
                action = :command
            end

        # Define the genes and choices for ArgParse
        genes = ["V", "D", "J"]
        choices = ["IGKV", "IGLV", "IGHV"]

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
            "-f", "--forwardarrayindex"
                help = "Forward array index to use for demultiplexing"
                arg_type = String
                default = ""
            "--case-filter-regex"
                help = "Regex pattern to filter cases (e.g., '[ACDERF]' for monkey cases)"
                arg_type = String
        end

        @add_arg_table! s["diff"] begin
            "fasta"
                help = "FASTA files with sequences"
                nargs = '+'
                required = true
        end

        @add_arg_table! s["linkage"] begin
        "input"
            help = "TSV/TSV.GZ with columns: case and db_name (or specify columns)"
            required = true
        "edges"
            help = "TSV.GZ to save pairwise edges table (metrics, p, q)"
            required = true
        "-C", "--case-col"
            help = "Name of column with donor/case id"
            default = "case"
            arg_type = String
        "-A", "--allele-col"
            help = "Name of column with allele name"
            default = "db_name"
            arg_type = String
        "-m", "--min-donors"
            help = "Minimum donors required to include an allele"
            default = 2
            arg_type = Int
            range_tester = (x->x >= 1)
        "--min-support"
            help = "Minimum co-present donors (n11) to consider an edge"
            default = 3
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-jaccard"
            help = "Minimum Jaccard co-presence required to consider an edge"
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--similarity"
            help = "Similarity mode: r or r2"
            default = "r"
            arg_type = String
            range_tester = (x->(x == "r") | (x == "r2"))
        "--threshold"
            help = "Similarity threshold for complete-linkage cut (on r or r2)"
            default = 0.5
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-cluster-size"
            help = "Minimum cluster size to output"
            default = 3
            arg_type = Int
            range_tester = (x->x >= 1)
        "--clusters"
            help = "Optional path to save clusters (TSV)"
            arg_type = String
        end

        @add_arg_table! s["haplotype"] begin
        "input"
            help = "TSV/TSV.GZ with columns: case and allele (or specify columns)"
            required = true
        "output"
            help = "TSV file to save haplotype inference results"
            required = true
        "-C", "--case-col"
            help = "Name of column with donor/case id"
            default = "case"
            arg_type = String
        "-A", "--allele-col"
            help = "Name of column with allele name"
            default = "db_name"
            arg_type = String
        "-G", "--gene-col"
            help = "Name of column with gene name (for grouping alleles by gene)"
            default = "gene"
            arg_type = String
        "-c", "--mincount"
            help = "Minimum count for an allele to be considered"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
        "-r", "--min-ratio"
            help = "Minimum allelic ratio threshold (minor/major allele ratio)"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-f", "--novel-fasta"
            help = "Optional FASTA file with novel alleles to mark in results"
            arg_type = String
        end

        @add_arg_table! s["table"] begin
            "outerjoin"
                help = "Perform outer join on two TSV files"
                action = :command
            "leftjoin"
                help = "Perform left join on two TSV files"
                action = :command
            "transform"
                help = "Transform TSV file using regex with capture groups"
                action = :command
            "aggregate"
                help = "Aggregate unique rows by selected columns"
                action = :command
            "unique"
                help = "Select unique rows based on specified columns"
                action = :command
            "sort"
                help = "Sort TSV file by specified columns"
                action = :command
            "filter"
                help = "Filter TSV file by column using regex or numeric operations"
                action = :command
            "select"
                help = "Select specific columns from TSV file"
                action = :command
        end

        @add_arg_table! s["table"]["outerjoin"] begin
            "left"
                help = "Left TSV file"
                required = true
                arg_type = String
            "right"
                help = "Right TSV file"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-k", "--keys"
                help = "Column names to join on (comma-separated)"
                required = true
                arg_type = String
            "--left-keys"
                help = "Column names to join on from left file (comma-separated, defaults to --keys)"
                arg_type = String
            "--right-keys"
                help = "Column names to join on from right file (comma-separated, defaults to --keys)"
                arg_type = String
            "--left-prefix"
                help = "Prefix for left file columns (defaults to no prefix)"
                arg_type = String
            "--right-prefix"
                help = "Prefix for right file columns (defaults to no prefix)"
                arg_type = String
            "--left-select"
                help = "Columns to select from left file (comma-separated, defaults to all columns)"
                arg_type = String
            "--right-select"
                help = "Columns to select from right file (comma-separated, defaults to all columns)"
                arg_type = String
        end

        @add_arg_table! s["table"]["leftjoin"] begin
            "left"
                help = "Left TSV file"
                required = true
                arg_type = String
            "right"
                help = "Right TSV file"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-k", "--keys"
                help = "Column names to join on (comma-separated)"
                required = true
                arg_type = String
            "--left-keys"
                help = "Column names to join on from left file (comma-separated, defaults to --keys)"
                arg_type = String
            "--right-keys"
                help = "Column names to join on from right file (comma-separated, defaults to --keys)"
                arg_type = String
            "--left-prefix"
                help = "Prefix for left file columns (defaults to no prefix)"
                arg_type = String
            "--right-prefix"
                help = "Prefix for right file columns (defaults to no prefix)"
                arg_type = String
            "--left-select"
                help = "Columns to select from left file (comma-separated, defaults to all columns)"
                arg_type = String
            "--right-select"
                help = "Columns to select from right file (comma-separated, defaults to all columns)"
                arg_type = String
        end

        @add_arg_table! s["table"]["transform"] begin
            "input"
                help = "Input TSV file to transform"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-c", "--column"
                help = "Column name to transform"
                required = true
                arg_type = String
            "-p", "--pattern"
                help = "Regex pattern with capture groups (e.g., 'ID_(\\d+)_(\\w+)')"
                required = true
                arg_type = String
            "-r", "--replacement"
                help = "Replacement string using capture groups (e.g., '\\1-\\2')"
                required = true
                arg_type = String
            "--new-column"
                help = "Optional name for new column containing captured groups"
                arg_type = String
        end

        @add_arg_table! s["table"]["aggregate"] begin
            "input"
                help = "Input TSV file to aggregate"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-g", "--group-by"
                help = "Column names to group by (comma-separated)"
                required = true
                arg_type = String
            "-k", "--keep-columns"
                help = "Additional columns to keep (comma-separated, defaults to all non-group columns)"
                arg_type = String
            "-c", "--count-column"
                help = "Name for the count column"
                default = "count"
                arg_type = String
        end

        @add_arg_table! s["table"]["unique"] begin
            "input"
                help = "Input TSV file to make unique"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Column names to select for uniqueness (comma-separated)"
                required = true
                arg_type = String
        end

        @add_arg_table! s["table"]["sort"] begin
            "input"
                help = "Input TSV file to sort"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Column names to sort by (comma-separated, in order of priority)"
                required = true
                arg_type = String
            "-r", "--reverse"
                help = "Sort in descending order (default: ascending)"
                action = :store_true
        end

        @add_arg_table! s["table"]["filter"] begin
            "input"
                help = "Input TSV file to filter"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-c", "--column"
                help = "Column name to filter on"
                required = true
                arg_type = String
            "--pattern"
                help = "Regex pattern for string filtering (use this for text columns)"
                arg_type = String
            "--operator"
                help = "Numeric operator: <, <=, >=, > (use with --threshold for numeric columns)"
                arg_type = String
                range_tester = (x->x ∈ ["<", "<=", ">=", ">"])
            "--threshold"
                help = "Numeric threshold value (use with --operator for numeric columns)"
                arg_type = Float64
        end

        @add_arg_table! s["table"]["select"] begin
            "input"
                help = "Input TSV file to select columns from"
                required = true
                arg_type = String
            "output"
                help = "Output TSV file"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Column names to select (comma-separated)"
                required = true
                arg_type = String
        end

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
                default = "IGHV"
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

        @add_arg_table! s["blast"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with database sequences"
            required = true
        "output"
            help = "TSV file to save discovery results"
            required = true
        "-p", "--pseudo"
            help = "FASTA file with pseudo-genes"
            arg_type = String
            default = ""
        "-c", "--minfullcount"
            help = "Minimum full cluster size"
            default = 5
            arg_type = Int
        "-f", "--minfullfreq"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-s", "--subjectcov"
            help = "Minimum subject (database) coverage"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-d", "--maxdist"
            help = "Maximum distance allowed for alleles"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 290
            arg_type = Int
            range_tester = (x->x >= 1)
        "-e", "--edge"
            help = "Minimum number of nucleotides required between target gene and end of the read"
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "-a", "--args"
            help = "Additional arguments to pass to blastn"
            arg_type = String
            default = "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
        "-o",  "--overwrite"
            help = "Overwrite existing files (i.e BLAST cache)"
            action = :store_true
        "-g", "--gene"
            help = "Use gene preset parameters for V, D, or J analysis (can be overridden by explicit parameters)"
            arg_type = String
            range_tester = (x->x ∈ keys(BLAST_PRESETS))
        "-G", "--show-presets"
            help = "Show preset parameters for V, D, and J analysis"
            action = :store_true
        "--forward"
            help = "Forward extension length"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "--reverse"
            help = "Reverse extension length"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "-q", "--minquality"
            help = "Minimum quality of the trimming alignment. Affixes with lower quality will drop the candidate."
            default = 0.75
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-i", "--isin"
            help = "Keep sequences that are substrings of known alleles"
            action = :store_false
        "--keep-failed"
            help = "Keep rows where trimming failed (aln_qseq empty). By default such rows are dropped."
            action = :store_true
        "--min-corecov"
            help = "Minimum ratio length(aln_qseq)/length(db_seq) after trimming"
            default = 0.6
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-v", "--verbose"
            help = "Print verbose output and save intermediate files"
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
        "-f", "--minratio"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c", "--mincount"
            help = "Minimum cluster size"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
        "--min-allele-mratio"
            help = "Minimum allelic ratio applied within each gene group for the allele against median"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "--min-gene-mratio"
            help = "Minimum allelic ratio applied within each gene group for the gene against median"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "--rss"
            help = "Comma-separated list of rss fragments: heptamer, spacer, nonamer"
            default = "heptamer"
            arg_type = String
        "--extension"
            help = "Length of extension on RSS side instead of RSS elements"
            arg_type = Int
        "--raw"
            help = "Unfiltered exact search results for diagnostics"
            arg_type = String
        "-n","--noplot"
            help = "Disable unicode gene plot"
            action = :store_true
        "-g", "--gene"
            default = "V"
            range_tester = (x->x ∈ genes)
            arg_type = String
            help = "gene; must be one of " * join(genes, ", ", " or ")
        "-a", "--affix"
            help = "Number of bases to extract from the non-RSS side of the sequence"
            arg_type = Int
            default = 13
            range_tester = (x->x >= 1)
        "-t", "--top"
            help = "Saves at most N records of flank and sequence."
            arg_type = Int
            default = 1
            range_tester = (x->x >= 1)
        "-r", "--refgene"
            help = "Space separated reference genes to use for computing ratio"
            nargs = '*'  # Accepts zero or more values
            arg_type = String
        "-l", "--limit"
            help = "Limit to this number of sequences, zero means no limit"
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "-e","--expect"
            help = "TSV file containing gene names and their corresponding allele_freq threshold, with two columns: name and ratio"
            arg_type = String
        "-d","--deletion"
            help = "TSV file containing gene names and their corresponding gene_count_freq threshold, with two columns: name and ratio"
            arg_type = String
        "--locus"
            help = "Locus to filter genes to start with this string (e.g. IGHV) excluding other genes from the analysis (i.e control genes)"
            arg_type = String
            default = "IG"
        "--ref-fasta"
            help = "Optional reference FASTA file to check if sequences are in database (adds isin_db column)"
            arg_type = String
        end

        @add_arg_table! s["nwpattern"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with aligned gene sequences to build trimming profile"
            required = true
        "output"
            help = "TSV file to save demultiplex data"
            required = true
        "-w", "--window"
            help = "Length of the flanks"
            default = 30
            arg_type = Int
            range_tester = (x->x > 1)
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 200
            arg_type = Int
            range_tester = (x->x >= 1)
        "-k", "--kmer"
            help = "Kmer size which will be used to search for patterns"
            default = 7
            arg_type = Int
            range_tester = (x->x >= 3)
        "-d", "--maxdist"
            help = "Maximum distance allowed for alleles"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 0)
        "-s", "--sample"
            help = "Number of kmers in a combination"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
        "--max-combinations"
            help = "Number of combinations to sample (higher slower but more sequences found)"
            default = 100000
            arg_type = Int
            range_tester = (x->x >= 1)
        "--max-attempts"
            help = "Number of samplings to attempt (higher slower but more sequences found)"
            default = 100000
            arg_type = Int
            range_tester = (x->x >= 1)
        "-f","--minfreq"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c","--mincount"
            help = "Minimum count for an allele"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
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
        "-t", "--top"
            help = "Save top candidates with highest counts per allele_name and case"
            arg_type = Int
            default = 5
            range_tester = (x->x >= 1)
        "-w", "--weights"
            help = "Length of the position weight matrix"
            default = 25
            arg_type = Int
            range_tester = (x->x > 1)
        "-n", "--noprofile"
            help = "Use distance to database gene lengths instead of profiles to trim reads"
            action = :store_true
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 200
            arg_type = Int
            range_tester = (x->x >= 1)
        "-k", "--kmer"
            help = "Kmer size which will be used to search for patterns"
            default = 12
            arg_type = Int
            range_tester = (x->x >= 3)
        "-u", "--usehalf"
            help = "Use half of the allele as a search word instead of kmer"
            action = :store_true
        "-m", "--maxkmer"
            help = "Maximum kmer size if kmer size needs to be increased automatically"
            default = 50
            arg_type = Int
            range_tester = (x->x >= 3)
        "-d", "--maxdist"
            help = "Maximum distance allowed for alleles"
            default = 50
            arg_type = Int
            range_tester = (x->x >= 0)
        "-s", "--sample"
            help = "Number of kmers to sample from the set of all kmers to search for a gene"
            default = 5
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

        @add_arg_table! s["regex"] begin
            "tsv"
                help = "TSV file with demultiplexed reads"
                required = true
            "fasta"
                help = "FASTA file with query alleles"
                required = true
            "output"
                help = "TSV file to save data with identified alleles"
                required = true
            "--insert-minlen"
                help = "Minimum length of an insert"
                default = 15
                arg_type = Int
                range_tester = (x->x >= 1)
            "--insert-maxlen"
                help = "Maximum length of an insert"
                default = 40
                arg_type = Int
                range_tester = (x->x >= 1)
            "--flank-mincount"
                help = "Minimum number of reads per flanks"
                default = 50
                arg_type = Int
                range_tester = (x->x >= 1)
            "--flank-frequency"
                help = "Lowest frequency of reads per flank"
                default = 0.5
                arg_type = Float64
                range_tester = (x-> (x >= 0.0) & (x <= 1.0))
            "--flank-minwell"
                help = "Minimum number of wells with reads per flank"
                default = 1
                arg_type = Int
                range_tester = (x->x >= 1)
            "-c", "--mincount"
                help = "Minimum count of a match"
                default = 3
                arg_type = Int
                range_tester = (x->x >= 1)
            "-d", "--maxdist"
                help = "Maximum distance of an insert to known allele"
                default = 3
                arg_type = Int
                range_tester = (x->x >= 1)
            "-f","--frequency"
                help = "Lowest frequency of a match"
                default = 0.1
                arg_type = Float64
                range_tester = (x-> (x >= 0.0) & (x <= 1.0))
            "-p", "--nprefix"
                help = "Number of nucleotides to extract from the 5' end of the query sequence"
                default = 7
                arg_type = Int
                range_tester = (x->x >= 1)
            "-s", "--nsuffix"
                help = "Number of nucleotides to extract from the 3' end of the query sequence"
                default = 7
                arg_type = Int
                range_tester = (x->x >= 1)
        end

        @add_arg_table! s["hsmm"] begin
        "tsv"
            help = "TSV/TSV.GZ demultiplex file with columns well, case, name, genomic_sequence"
            required = true
        "fasta"
            help = "FASTA file with D alleles (known reference)"
            required = true
        "output"
            help = "TSV.GZ file to save detected D alleles (novel and/or known) with flanks"
            required = true
        "-r", "--ratio"
            help = "Allelic ratio threshold for known D selection per donor and gene"
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c", "--mincount"
            help = "Minimum count for a known D allele to be considered in training"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "--min-posterior"
            help = "Minimum posterior probability for accepting an HSMM detection"
            default = 0.7
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-gene-len"
            help = "Minimum D gene length for HSMM duration model (auto if 0)"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 0)
        "--max-gene-len"
            help = "Maximum D gene length for HSMM duration model (auto if 0)"
            default = 70
            arg_type = Int
            range_tester = (x->x >= 0)
        "-l", "--limit"
            help = "Limit number of demultiplexed reads to process (0 means no limit)"
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--out-mincount"
            help = "Minimum count for an extracted D (after HSMM) to keep in output"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "--out-minratio"
            help = "Minimum allelic ratio within gene (per donor) for output filtering"
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-heptamer-prob-pre"
            help = "Minimum probability under pre-heptamer PWM to keep detection (0 disables)"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-heptamer-prob-post"
            help = "Minimum probability under post-heptamer PWM to keep detection (0 disables)"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        @add_arg_table! s["fasta"] begin
        "input"
            help = "TSV file with data to extract"
            required = true
        "output"
            help = "Output FASTA file"
            required = true
        "-n", "--colname"
            help = "Name of the column with sequence names/IDs"
            default = "allele_name"
            arg_type = String
        "-s", "--colseq"
            help = "Name of the column with sequences"
            default = "seq"
            arg_type = String
        "-d", "--coldesc"
            help = "Optional name of the column with descriptions (used in FASTA headers)"
            default = nothing
            arg_type = Union{String, Nothing}
        "-f", "--filter"
            help = "Optional regex pattern to filter colname column (e.g., 'Novel')"
            default = nothing
            arg_type = Union{String, Nothing}
        "-c", "--cleanup"
            help = "Optional regex pattern to remove from sequence names (e.g., ' Novel')"
            default = nothing
            arg_type = Union{String, Nothing}
        "--desc-filter"
            help = "Optional regex pattern to filter coldesc column (without capture groups: filters only; with capture groups: filters and includes captured string in FASTA header, e.g., '(case\\d+)')"
            default = nothing
            arg_type = Union{String, Nothing}
        "--no-sort"
            help = "Disable sorting alleles by name (default: sort enabled)"
            action = :store_true
        "--mincase"
            help = "Minimum number of donors (cases) that must have this allele to include it in output"
            default = 1
            arg_type = Int
            range_tester = (x->x >= 1)
        "--case-col"
            help = "Name of the column containing donor/case identifiers"
            default = "case"
            arg_type = String
        end

        @add_arg_table! s["merge"] begin
        "output"
            help = "Output merged FASTA file"
            required = true
        "inputs"
            help = "Input FASTA files to merge (2 or more files)"
            nargs = '+'
            required = true
        "-c", "--cleanup"
            help = "Optional regex pattern to remove from sequence names (e.g., ' Novel')"
            default = nothing
            arg_type = Union{String, Nothing}
        "--no-sort"
            help = "Disable sorting sequences by name (default: sort enabled)"
            action = :store_true
        "--prefer-last"
            help = "When duplicate sequences have different names, prefer the last encountered (default: prefer first)"
            action = :store_true
        "--add-source-prefix"
            help = "Add source filename as prefix to sequence names"
            action = :store_true
        end

        @add_arg_table! s["collect"] begin
        "pattern"
            help = "Pattern of TSV files to collect"
            required = true
        "output"
            help = "TSV file to save data"
            required = true
        end

        @add_arg_table! s["hash"] begin
        "fastain"
            help = "FASTA file with allele names"
            required = true
        end

        @add_arg_table! s["exclude"] begin
        "input"
            help = "TSV file with columns allele_name and seq"
            required = true
        "output"
            help = "TSV file to save filtered input"
            required = true
        "fasta"
            help = "FASTA file with allele names and sequences to exclude"
            required = true
        "-n", "--colname"
            help = "Name of the column with allele names"
            default = "allele_name"
            arg_type = String
        "-s", "--colseq"
            help = "Name of the column with sequence"
            default = "seq"
            arg_type = String
        end

        @add_arg_table! s["bwa"] begin
        "tsv"
            help = "TSV file with columns allele_name and seq"
            required = true
        "output"
            help = "TSV file to save filtered input"
            required = true
        "genome"
            help = "FASTA file with indexed genome"
            required = true
            nargs='+'
            arg_type = String
        "-c", "--chromosome"
            help = "Chromosome string to filter by"
            default = "chromosome 14"
            arg_type = String
        "-n", "--colname"
            help = "Name of the column with allele names"
            default = "best_name"
            arg_type = String
        "-s", "--colseq"
            help = "List column names with sequences"
            default = ["prefix", "best_aln", "suffix"]
            nargs = '*'  # Accepts zero or more values
            arg_type = String
        "-t", "--tag"
            help = "Regex to filter valid descriptions of chromosomes"
            default = "(.*Primary Assembly.*)|(.*alternate locus.*)"
            arg_type = String
        end

        # Show help if no arguments are provided
        try
            open("immunediscover.log", "a") do io
                logger = ConsoleLogger(io)
                with_logger(logger) do
                    @info "$SOFTWARE_VERSION $(Dates.now()) - Parsing command line arguments: $args"
                end
            end
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
