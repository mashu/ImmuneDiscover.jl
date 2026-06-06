function add_table_args!(s)
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

        for (cmd, tag) in (("outerjoin", "outerjoin"), ("leftjoin", "leftjoin"))
            jn = s["table"][cmd]
            @add_arg_table! jn begin
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
            end

            add_arg_group!(jn, "Join keys", tag * "_keys")
            @add_arg_table! jn begin
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
            end

            add_arg_group!(jn, "Column selection", tag * "_cols")
            @add_arg_table! jn begin
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

        fl = s["table"]["filter"]
        @add_arg_table! fl begin
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
        end

        add_arg_group!(fl, "String filter", "filter_string")
        @add_arg_table! fl begin
            "--pattern"
                help = "Regex pattern for string filtering (for text columns)"
                arg_type = String
        end

        add_arg_group!(fl, "Numeric filter", "filter_numeric")
        @add_arg_table! fl begin
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

        tf = s["table"]["fasta"]
        @add_arg_table! tf begin
        "input"
            help = "Input TSV file path"
            required = true
        "output"
            help = "Output FASTA file path"
            required = true
        end

        add_arg_group!(tf, "Input columns", "tfasta_cols")
        @add_arg_table! tf begin
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
        end

        add_arg_group!(tf, "Name and description filtering", "tfasta_filter")
        @add_arg_table! tf begin
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
        end

        add_arg_group!(tf, "Case filtering", "tfasta_case")
        @add_arg_table! tf begin
        "--mincase"
            help = "Minimum number of cases that must include the allele to export"
            default = 1
            arg_type = Int
            range_tester = (x->x >= 1)
        "--case-col"
            help = "Column name with donor/case identifiers"
            default = "case"
            arg_type = String
        end

        add_arg_group!(tf, "Output", "tfasta_out")
        @add_arg_table! tf begin
        "--no-sort"
            help = "Do not sort records by sequence name"
            action = :store_true
        "--unique-sequences"
            help = "Keep only unique sequences (ignore sequence names, use first encountered name)"
            action = :store_true
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
    return s
end
