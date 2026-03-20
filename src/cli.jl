module Cli
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
        "minfullratio" => 0.1,
        "length" => 290,
        "maxdist" => 20,
        "minfullcount" => 5,
        "edge" => 0,
        "subjectcov" => 0.1,
        "minquality" => 0.75,
        "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
    )

    # V preset tuned for genomic IGHV novel recovery (see docs + tuning/README): more BLAST hits per
    # read, slightly relaxed trim/core coverage and allelic ratio, longer extensions for affix alignment.
    const BLAST_PRESETS = Dict(
        "V" => Dict(
            "forward" => 20,
            "reverse" => 20,
            "minfullratio" => 0.035,
            "length" => 283,
            "maxdist" => 14,
            "minfullcount" => 5,
            "minquality" => 0.62,
            "min-corecov" => 0.50,
            "args" => "-task megablast -subject_besthit -num_alignments 25 -qcov_hsp_perc 50"
        ),
        "D" => Dict(
            "forward" => 40,
            "reverse" => 40,
            "minfullratio" => 0.2,
            "length" => 5,
            "maxdist" => 20,
            "minfullcount" => 10,
            "edge" => 10,
            "subjectcov" => 0.25,
            "minquality" => 0.5,
            "args" => "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5"
        ),
        "J" => Dict(
            "forward" => 12,
            "reverse" => 12,
            "minfullratio" => 0.1,
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

    function get_blast_block(args)
        cmd = get(args, "%COMMAND%", "")
        cmd == "blast" && return args["blast"]
        cmd == "discover" && get(args["discover"], "%COMMAND%", "") == "blast" && return args["discover"]["blast"]
        return nothing
    end

    function show_blast_params(args)
        block = get_blast_block(args)
        block === nothing && return
        gene = block["gene"]
        if haskey(BLAST_PRESETS, gene)
            preset = BLAST_PRESETS[gene]
            println("BLAST Preset for $gene gene:")
            for (param, value) in preset
                println("  --$param = $(block[param])")
            end
        end
    end

    function apply_blast_presets!(parsed_args)
        block = get_blast_block(parsed_args)
        block === nothing && return parsed_args
        gene = block["gene"]
        if haskey(BLAST_PRESETS, gene)
            preset = BLAST_PRESETS[gene]

            # Only apply preset if the current value equals the CLI default
            # This means the user didn't explicitly override it
            for (key, preset_value) in preset
                if haskey(block, key) && haskey(BLAST_CLI_DEFAULTS, key)
                    current_value = block[key]
                    default_value = BLAST_CLI_DEFAULTS[key]
                    
                    # Apply preset only if current value equals the default
                    if current_value == default_value
                        @info "Applying $gene preset $key: $current_value → $preset_value"
                        block[key] = preset_value
                    else
                        @info "Keeping user override for $key: $current_value (not default $default_value)"
                    end
                else
                    # Apply preset for keys not in CLI defaults
                    if haskey(block, key)
                        @info "Applying $gene preset $key: $(get(block, key, nothing)) → $preset_value"
                        block[key] = preset_value
                    end
                end
            end
        end
        return parsed_args
    end

    # Lazy version detection — avoids running external `git` at const-initialization
    # time, which triggers a Julia 1.12 compiler inference bug.
    const version_cache = Ref{String}("")
    const hash_cache = Ref{String}("")

    function run_git_or_unknown(args::Cmd)
        Sys.which("git") === nothing && return "unknown"
        proc = run(pipeline(args, stderr=devnull), wait=false)
        wait(proc)
        success(proc) || return "unknown"
        return strip(read(args, String))
    end

    function software_version()
        if isempty(version_cache[])
            version_cache[] = run_git_or_unknown(`git -C $(@__DIR__) describe --tags --abbrev=0`)
        end
        return version_cache[]
    end

    function software_git_hash()
        if isempty(hash_cache[])
            hash_cache[] = run_git_or_unknown(`git -C $(@__DIR__) rev-parse HEAD`)
        end
        return hash_cache[]
    end

    """
        always_gz(file_path)

    Return path that ends with .gz, appending the extension if needed.
    """
    function always_gz(file_path)
        endswith(file_path, ".gz") ? file_path : file_path * ".gz"
    end

    """
        parse_commandline(args)

    Handle command line
    """
    function parse_commandline(args)
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "$(software_version()) (git $(software_git_hash()))",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "GKHLab, $(software_version()) (git $(software_git_hash()))")
        @add_arg_table! s begin
            "discover"
                help = "De novo allele discovery (blast, hsmm)"
                action = :command
            "search"
                help = "Search against known references (exact, heptamer, bwa)"
                action = :command
            "analyze"
                help = "Downstream analysis (cooccurrence, haplotype)"
                action = :command
            "preprocess"
                help = "Data preparation (demultiplex)"
                action = :command
            "table"
                help = "Table utilities (outerjoin, leftjoin, transform, aggregate, unique, sort, filter, select, fasta, collect, exclude)"
                action = :command
            "fasta"
                help = "FASTA utilities (merge, diff, hash)"
                action = :command
            end

        # Define the genes and choices for ArgParse
        genes = ["V", "D", "J"]
        choices = ["IGKV", "IGLV", "IGHV"]

        @add_arg_table! s["preprocess"] begin
            "demultiplex"
                help = "Demultiplex indexed plate libraries into a TSV with per-read metadata"
                action = :command
        end

        @add_arg_table! s["preprocess"]["demultiplex"] begin
            "fastq"
                help = "Input FASTQ file with reads (single-end)"
                required = true
            "indices"
                help = "TSV with well indices/barcodes and case identifiers"
                required = true
            "output"
                help = "Output TSV (gz auto-enabled) with demultiplexed reads and metadata"
                required = true
            "-l", "--length"
                help = "Minimum read length to keep"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 200
            "-s", "--split"
                help = "Write per-case FASTQ files"
                action = :store_true
            "-f", "--forwardarrayindex"
                help = "Name of the forward array index to use for demultiplexing (if present)"
                arg_type = String
                default = ""
            "--case-filter-regex"
                help = "Regex to keep only cases matching pattern (e.g., '[ACDERF]')"
                arg_type = String
        end

        @add_arg_table! s["analyze"] begin
            "cooccurrence"
                help = "Analyze cross-donor co-occurrence among alleles across cases"
                action = :command
            "haplotype"
                help = "Infer approximate haplotypes per case using diploid assumptions"
                action = :command
        end

        # Options for cooccurrence (canonical)
        @add_arg_table! s["analyze"]["cooccurrence"] begin
        "input"
            help = "TSV/TSV.GZ with columns: case and db_name (or specify columns)"
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
        "--cluster-method"
            help = "Clustering method on rho: components, complete, average, or single"
            default = "components"
            arg_type = String
            range_tester = (x-> (x ∈ ["components","complete","average","single"]))
        "--cluster-threshold"
            help = "Similarity threshold for complete-linkage clustering on rho (0..1)"
            default = 0.5
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--debug-triangles"
            help = "Print rho-based triangle diagnostics at cluster-threshold"
            action = :store_true
        "--min-cluster-size"
            help = "Minimum cluster size to output"
            default = 3
            arg_type = Int
            range_tester = (x->x >= 1)
        "--clusters"
            help = "Optional path to save clusters (TSV)"
            arg_type = String
        end

        @add_arg_table! s["analyze"]["haplotype"] begin
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
                help = "Outer-join two TSVs on specified key columns"
                action = :command
            "leftjoin"
                help = "Left-join two TSVs on specified key columns"
                action = :command
            "transform"
                help = "Regex-based column transform with capture-group replacement and optional new column"
                action = :command
            "aggregate"
                help = "Group by selected columns and count unique groups (keep optional columns)"
                action = :command
            "unique"
                help = "Select distinct rows by specified columns"
                action = :command
            "sort"
                help = "Sort TSV by one or more columns (asc/desc)"
                action = :command
            "filter"
                help = "Filter rows by string regex or numeric threshold operation"
                action = :command
            "select"
                help = "Project a subset of columns from a TSV"
                action = :command
            "fasta"
                help = "Export sequences from TSV to FASTA with optional filtering/cleanup"
                action = :command
            "collect"
                help = "Concatenate many TSVs into a single table (schema must match)"
                action = :command
            "exclude"
                help = "Exclude sequences whose names or sequences overlap with a FASTA reference"
                action = :command
        end

        @add_arg_table! s["table"]["outerjoin"] begin
            "left"
                help = "Left TSV file path"
                required = true
                arg_type = String
            "right"
                help = "Right TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-k", "--keys"
                help = "Comma-separated column names to join on"
                required = true
                arg_type = String
            "--left-keys"
                help = "Comma-separated join keys from left file (defaults to --keys)"
                arg_type = String
            "--right-keys"
                help = "Comma-separated join keys from right file (defaults to --keys)"
                arg_type = String
            "--left-prefix"
                help = "Optional prefix for left non-key columns"
                arg_type = String
            "--right-prefix"
                help = "Optional prefix for right non-key columns"
                arg_type = String
            "--left-select"
                help = "Comma-separated subset of columns to keep from left file"
                arg_type = String
            "--right-select"
                help = "Comma-separated subset of columns to keep from right file"
                arg_type = String
        end

        @add_arg_table! s["table"]["leftjoin"] begin
            "left"
                help = "Left TSV file path"
                required = true
                arg_type = String
            "right"
                help = "Right TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-k", "--keys"
                help = "Comma-separated column names to join on"
                required = true
                arg_type = String
            "--left-keys"
                help = "Comma-separated join keys from left file (defaults to --keys)"
                arg_type = String
            "--right-keys"
                help = "Comma-separated join keys from right file (defaults to --keys)"
                arg_type = String
            "--left-prefix"
                help = "Optional prefix for left non-key columns"
                arg_type = String
            "--right-prefix"
                help = "Optional prefix for right non-key columns"
                arg_type = String
            "--left-select"
                help = "Comma-separated subset of columns to keep from left file"
                arg_type = String
            "--right-select"
                help = "Comma-separated subset of columns to keep from right file"
                arg_type = String
        end

        @add_arg_table! s["table"]["transform"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-c", "--column"
                help = "Target column to transform"
                required = true
                arg_type = String
            "-p", "--pattern"
                help = "Regex with capture groups (e.g., 'ID_(\\d+)_(\\w+)')"
                required = true
                arg_type = String
            "-r", "--replacement"
                help = "Replacement string using capture groups (e.g., '\\1-\\2')"
                required = true
                arg_type = String
            "--new-column"
                help = "Optional name for a new column storing captured groups"
                arg_type = String
        end

        @add_arg_table! s["table"]["aggregate"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-g", "--group-by"
                help = "Comma-separated column names to group by"
                required = true
                arg_type = String
            "-k", "--keep-columns"
                help = "Comma-separated additional columns to keep (defaults to all non-group columns)"
                arg_type = String
            "-c", "--count-column"
                help = "Name of the count column"
                default = "count"
                arg_type = String
        end

        @add_arg_table! s["table"]["unique"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Comma-separated columns to form distinct rows"
                required = true
                arg_type = String
        end

        @add_arg_table! s["table"]["sort"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Comma-separated columns to sort by (priority order)"
                required = true
                arg_type = String
            "-r", "--reverse"
                help = "Sort in descending order (default ascending)"
                action = :store_true
        end

        @add_arg_table! s["table"]["filter"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-c", "--column"
                help = "Target column to filter on"
                required = true
                arg_type = String
            "--pattern"
                help = "Regex pattern for string filtering (for text columns)"
                arg_type = String
            "--operator"
                help = "Numeric operator: <, <=, >=, > (requires --threshold)"
                arg_type = String
                range_tester = (x->x ∈ ["<", "<=", ">=", ">"])
            "--threshold"
                help = "Numeric threshold value (used with --operator for numeric columns)"
                arg_type = Float64
        end

        @add_arg_table! s["table"]["select"] begin
            "input"
                help = "Input TSV file path"
                required = true
                arg_type = String
            "output"
                help = "Output TSV (gz auto-enabled)"
                required = true
                arg_type = String
            "-c", "--columns"
                help = "Comma-separated list of columns to select"
                required = true
                arg_type = String
        end

        @add_arg_table! s["discover"] begin
            "blast"
                help = "BLAST-based candidate discovery with trimming, filtering, and identity clustering"
                action = :command
            "hsmm"
                help = "Detect D genes using an HSMM trained on RSS flanks (V/J masked)"
                action = :command
        end

        @add_arg_table! s["search"] begin
            "exact"
                help = "Exact match search of reads to database alleles with robust filters"
                action = :command
            "heptamer"
                help = "Identify heptamer RSS positions and extend/trim V reads accordingly"
                action = :command
            "bwa"
                help = "Genome mapping QC with BWA to retain sequences mapped to the target chromosome"
                action = :command
        end

        @add_arg_table! s["search"]["bwa"] begin
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

        @add_arg_table! s["search"]["heptamer"] begin
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

        @add_arg_table! s["discover"]["blast"] begin
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
        "-f", "--minfullratio"
            help = "Minimum allelic ratio within each gene group (count / max in gene)"
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
            help = "Minimum fraction (0–1) of affix positions that match the read in the semi-global affix–read alignment used for 5'/3' trimming; prefix and suffix each must meet this or the candidate is dropped."
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

        @add_arg_table! s["search"]["exact"] begin
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
        "--border"
            help = "Number of nucleotides from both read ends forming a border; reject if extension overlaps"
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "--adjust-per-gene-extension"
            help = "Auto-reduce extension per gene and side to avoid crossing the border"
            action = :store_true
        "--adjust-percent"
            help = "Target fraction (0-1] of reads per gene that should be safe (no border overlap) when calibrating per-gene extension"
            arg_type = Float64
            default = 1.0
            range_tester = (x-> (x > 0.0) & (x <= 1.0))
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
            help = "TSV file containing gene names and their corresponding gene_case_freq threshold, with two columns: name and ratio"
            arg_type = String
        "--locus"
            help = "Locus to filter genes to start with this string (e.g. IGHV) excluding other genes from the analysis (i.e control genes)"
            arg_type = String
            default = "IG"
        "--ref-fasta"
            help = "Optional reference FASTA file to check if sequences are in database (adds isin_db column)"
            arg_type = String
        end

        @add_arg_table! s["discover"]["hsmm"] begin
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

        # Table → fasta export
        @add_arg_table! s["table"]["fasta"] begin
        "input"
            help = "Input TSV file path"
            required = true
        "output"
            help = "Output FASTA file path"
            required = true
        "-n", "--colname"
            help = "Column name with sequence names/IDs"
            default = "allele_name"
            arg_type = String
        "-s", "--colseq"
            help = "Column name(s) with nucleotide sequences (comma-separated for concatenation)"
            default = "seq"
            arg_type = String
        "-d", "--coldesc"
            help = "Optional column with descriptions appended to FASTA headers"
            default = nothing
            arg_type = Union{String, Nothing}
        "-f", "--filter"
            help = "Regex to filter sequence names (e.g., 'Novel')"
            default = nothing
            arg_type = Union{String, Nothing}
        "-c", "--cleanup"
            help = "Regex to remove from sequence names (e.g., ' Novel')"
            default = nothing
            arg_type = Union{String, Nothing}
        "--desc-filter"
            help = "Regex to filter description column; capture group 1 (if present) is appended"
            default = nothing
            arg_type = Union{String, Nothing}
        "--no-sort"
            help = "Do not sort records by sequence name"
            action = :store_true
        "--mincase"
            help = "Minimum number of cases that must include the allele to export"
            default = 1
            arg_type = Int
            range_tester = (x->x >= 1)
        "--case-col"
            help = "Column name with donor/case identifiers"
            default = "case"
            arg_type = String
        "--unique-sequences"
            help = "Keep only unique sequences (ignore sequence names, use first encountered name)"
            action = :store_true
        end

        # FASTA group (operations on FASTA files)
        @add_arg_table! s["fasta"] begin
            "merge"
                help = "Merge multiple FASTA files, keeping only unique sequences"
                action = :command
            "diff"
                help = "Diff two FASTA files based on sequence identity but keep associated names"
                action = :command
            "hash"
                help = "Add hash based _S suffix to all allele names in the FASTA file"
                action = :command
        end

        @add_arg_table! s["fasta"]["merge"] begin
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

        @add_arg_table! s["fasta"]["diff"] begin
            "fasta"
                help = "FASTA files with sequences"
                nargs = '+'
                required = true
        end

        @add_arg_table! s["fasta"]["hash"] begin
        "fastain"
            help = "Input FASTA file path"
            required = true
        end

        @add_arg_table! s["table"]["collect"] begin
        "pattern"
            help = "Glob pattern of TSV files to concatenate (columns must match)"
            required = true
        "output"
            help = "Output TSV file path for concatenated table"
            required = true
        end

        @add_arg_table! s["table"]["exclude"] begin
        "input"
            help = "Input TSV with allele_name and seq columns"
            required = true
        "output"
            help = "Output TSV with excluded rows removed"
            required = true
        "fasta"
            help = "Reference FASTA with names/sequences to exclude against"
            required = true
        "-n", "--colname"
            help = "Column name containing allele names"
            default = "allele_name"
            arg_type = String
        "-s", "--colseq"
            help = "Column name containing nucleotide sequences"
            default = "seq"
            arg_type = String
        end

        # Log command invocation (non-critical, silently skip on any failure)
        let logpath = "immunediscover.log"
            io = nothing
            opened = false
            if isdir(dirname(abspath(logpath)))
                io = open(logpath, "a")
                opened = true
            end
            if opened
                logger = ConsoleLogger(io)
                with_logger(logger) do
                    @info "$(software_version()) $(Dates.now()) - Parsing command line arguments: $args"
                end
                close(io)
            end
        end

        # CLI-boundary catch: ArgParse throws ArgParseError by design for invalid user input.
        # This is the standard pattern and acceptable at the CLI entry point.
        try
            return parse_args(args, s)
        catch e
            if e isa ArgParseError
                println(e)
                ArgParse.show_help(s)
            else
                rethrow()
            end
        end
    end
end
