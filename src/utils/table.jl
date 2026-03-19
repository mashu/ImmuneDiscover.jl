module Table
    using CSV
    using DataFrames
    using Logging
    using Glob

    export outerjoin_tsv, leftjoin_tsv, filter_tsv, transform_tsv, aggregate_tsv, unique_tsv, sort_tsv, select_tsv
    export handle_table

    # --- Unified join implementation ---

    """
        join_tsv(left_file, right_file, output_file; how, left_keys, right_keys, left_prefix, right_prefix, left_select, right_select)

    Perform a join (outer or left) on two TSV files by specified key columns.
    """
    function join_tsv(left_file, right_file, output_file;
                      how::Symbol,
                      left_keys, right_keys=nothing, left_prefix=nothing, right_prefix=nothing,
                      left_select=nothing, right_select=nothing)

        @info "Loading left file: $left_file"
        left_df = CSV.File(left_file, delim='\t') |> DataFrame

        @info "Loading right file: $right_file"
        right_df = CSV.File(right_file, delim='\t') |> DataFrame

        if right_keys === nothing
            right_keys = left_keys
        end
        if left_prefix === nothing
            left_prefix = ""
        end
        if right_prefix === nothing
            right_prefix = ""
        end

        if left_select !== nothing
            @info "Selecting columns from left file: $left_select"
            left_df = left_df[:, left_select]
        end
        if right_select !== nothing
            @info "Selecting columns from right file: $right_select"
            right_df = right_df[:, right_select]
        end

        @info "Left file: $(nrow(left_df)) rows, $(ncol(left_df)) columns"
        @info "Right file: $(nrow(right_df)) rows, $(ncol(right_df)) columns"
        @info "Join keys - Left: $left_keys, Right: $right_keys"

        left_df_prefixed = copy(left_df)
        right_df_prefixed = copy(right_df)

        for col in names(left_df)
            if !(col in left_keys) && !isempty(left_prefix)
                rename!(left_df_prefixed, col => "$(left_prefix)_$(col)")
            end
        end
        for col in names(right_df)
            if !(col in right_keys) && !isempty(right_prefix)
                rename!(right_df_prefixed, col => "$(right_prefix)_$(col)")
            end
        end

        JOIN_FNS = Dict{Symbol, Function}(:outer => outerjoin, :left => leftjoin)
        join_fn = get(JOIN_FNS, how) do
            error("Unsupported join type: $how. Use :outer or :left.")
        end
        @info "Performing $how join"
        result_df = join_fn(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)

        @info "Result: $(nrow(result_df)) rows, $(ncol(result_df)) columns"

        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, result_df, compress=true, delim='\t')
        @info "Output saved to: $output_gz"
        return result_df
    end

    function outerjoin_tsv(left_file, right_file, output_file; kwargs...)
        join_tsv(left_file, right_file, output_file; how=:outer, kwargs...)
    end

    function leftjoin_tsv(left_file, right_file, output_file; kwargs...)
        join_tsv(left_file, right_file, output_file; how=:left, kwargs...)
    end

    is_numeric_value(::Missing) = false
    is_numeric_value(::Number) = true
    is_numeric_value(::Any) = false

    function is_numeric_column(col)
        for val in col
            ismissing(val) && continue
            return is_numeric_value(val)
        end
        return false
    end

    to_float64(::Missing) = nothing
    to_float64(x::Number) = Float64(x)
    to_float64(x) = tryparse(Float64, string(x))

    parse_numeric_value(val) = to_float64(val)

    """
        filter_tsv(input_file, output_file, column; pattern=nothing, operator=nothing, threshold=nothing)

    Filter TSV file by column using regex or numeric operations.
    """
    function filter_tsv(input_file, output_file, column; pattern=nothing, operator=nothing, threshold=nothing)

        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame

        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"

        if !(column in names(df))
            error("Column '$column' not found in input file. Available columns: $(join(names(df), ", "))")
        end

        if pattern !== nothing && (operator !== nothing || threshold !== nothing)
            error("Cannot use both regex pattern and numeric operator/threshold. Choose one filtering method.")
        elseif pattern !== nothing
            @info "Filtering column '$column' using regex pattern: $pattern"
            regex = Regex(pattern)
            filter_mask = map(x -> ismissing(x) ? false : occursin(regex, string(x)), df[!, column])
        elseif operator !== nothing && threshold !== nothing
            col_data = df[!, column]
            if is_numeric_column(col_data)
                numeric_values = col_data
            else
                parsed = [parse_numeric_value(x) for x in col_data]
                if any(isnothing, parsed)
                    error("Column '$column' cannot be converted to numeric values for filtering")
                end
                numeric_values = [p for p in parsed]
            end

            @info "Filtering column '$column' using numeric operation: $operator $threshold"

            COMPARE_OPS = Dict{String, Function}(
                "<"  => (x, t) -> ismissing(x) ? false : x < t,
                "<=" => (x, t) -> ismissing(x) ? false : x <= t,
                ">=" => (x, t) -> ismissing(x) ? false : x >= t,
                ">"  => (x, t) -> ismissing(x) ? false : x > t,
            )
            compare_fn = get(COMPARE_OPS, operator) do
                error("Invalid operator: $operator. Must be one of: <, <=, >=, >")
            end
            filter_mask = map(x -> compare_fn(x, threshold), numeric_values)
        else
            error("Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering")
        end

        filtered_df = df[filter_mask, :]
        @info "Filter results: $(nrow(filtered_df)) rows kept, $(nrow(df) - nrow(filtered_df)) rows removed"

        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, filtered_df, compress=true, delim='\t')
        @info "Filtered data saved to: $output_gz"
        return filtered_df
    end

    function transform_tsv(input_file, output_file; column, pattern, replacement, new_column=nothing)
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"

        if !(column in names(df))
            error("Column '$column' not found in input file. Available columns: $(join(names(df), ", "))")
        end

        @info "Transforming column '$column' using pattern: $pattern"
        regex = Regex(pattern)
        original_values = df[!, column]
        transformed_values = String[]
        extracted_values = String[]

        for value in original_values
            m = match(regex, string(value))
            if m !== nothing
                result = replacement
                for i in 1:length(m.captures)
                    if m.captures[i] !== nothing
                        result = replace(result, "\\$i" => m.captures[i])
                    end
                end
                push!(transformed_values, result)
                if new_column !== nothing
                    push!(extracted_values, join([c for c in m.captures if c !== nothing], ""))
                else
                    push!(extracted_values, "")
                end
            else
                push!(transformed_values, string(value))
                push!(extracted_values, "")
            end
        end

        changed_count = sum(string.(original_values) .!= transformed_values)
        @info "Transformed $changed_count values in column '$column'"
        df[!, column] = transformed_values
        if new_column !== nothing
            df[!, new_column] = extracted_values
            @info "Added new column '$new_column'"
        end

        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, df, compress=true, delim='\t')
        @info "Transformed data saved to: $output_gz"
        return df
    end

    function aggregate_tsv(input_file, output_file; group_by, keep_columns=nothing, count_column="count")
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"

        for col in group_by
            if !(col in names(df))
                error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
            end
        end

        if keep_columns === nothing
            keep_columns = [col for col in names(df) if !(col in group_by)]
        else
            for col in keep_columns
                if !(col in names(df))
                    error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
                end
            end
        end

        selected_columns = vcat(group_by, keep_columns)
        df_selected = df[:, selected_columns]
        grouped_df = combine(groupby(df_selected, group_by), nrow => count_column)

        if !isempty(keep_columns)
            for col in keep_columns
                if !(col in names(grouped_df))
                    first_values = combine(groupby(df_selected, group_by), col => first => col)
                    grouped_df = leftjoin(grouped_df, first_values, on=group_by)
                end
            end
        end

        @info "Aggregated from $(nrow(df)) rows to $(nrow(grouped_df)) unique groups"
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, grouped_df, compress=true, delim='\t')
        return grouped_df
    end

    function unique_tsv(input_file, output_file; columns)
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found. Available: $(join(names(df), ", "))")
            end
        end
        unique_df = unique(df[:, columns])
        @info "Unique rows: $(nrow(unique_df)) (removed $(nrow(df) - nrow(unique_df)) duplicates)"
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, unique_df, compress=true, delim='\t')
        return unique_df
    end

    function sort_tsv(input_file, output_file; columns, reverse=false)
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found. Available: $(join(names(df), ", "))")
            end
        end
        sort!(df, columns, rev=reverse)
        @info "Sorted $(nrow(df)) rows by $columns ($(reverse ? "desc" : "asc"))"
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, df, compress=true, delim='\t')
        return df
    end

    function select_tsv(input_file, output_file; columns)
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found. Available: $(join(names(df), ", "))")
            end
        end
        selected_df = df[:, columns]
        @info "Selected $(ncol(selected_df)) columns from $(ncol(df))"
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, selected_df, compress=true, delim='\t')
        return selected_df
    end

    # --- Subcommand dispatch table for table operations ---

    function route_join(parsed_args, subcmd, how)
        block = parsed_args["table"][subcmd]
        left_keys = split(block["keys"], ",")
        right_keys = get(block, "right-keys", nothing)
        if right_keys !== nothing; right_keys = split(right_keys, ","); end
        left_prefix = get(block, "left-prefix", nothing)
        right_prefix = get(block, "right-prefix", nothing)
        left_select = get(block, "left-select", nothing)
        right_select = get(block, "right-select", nothing)
        if left_select !== nothing; left_select = split(left_select, ","); end
        if right_select !== nothing; right_select = split(right_select, ","); end

        join_tsv(block["left"], block["right"], block["output"];
                 how=how, left_keys=left_keys, right_keys=right_keys,
                 left_prefix=left_prefix, right_prefix=right_prefix,
                 left_select=left_select, right_select=right_select)
    end

    function handle_outerjoin(parsed_args, _, _)
        @info "Performing outer join"
        route_join(parsed_args, "outerjoin", :outer)
    end

    function handle_leftjoin(parsed_args, _, _)
        @info "Performing left join"
        route_join(parsed_args, "leftjoin", :left)
    end

    function handle_transform(parsed_args, _, _)
        @info "Transforming TSV file"
        new_column = get(parsed_args["table"]["transform"], "new-column", nothing)
        transform_tsv(
            parsed_args["table"]["transform"]["input"],
            parsed_args["table"]["transform"]["output"];
            column=parsed_args["table"]["transform"]["column"],
            pattern=parsed_args["table"]["transform"]["pattern"],
            replacement=parsed_args["table"]["transform"]["replacement"],
            new_column=new_column
        )
    end

    function handle_aggregate(parsed_args, _, _)
        @info "Aggregating TSV file"
        group_by = split(parsed_args["table"]["aggregate"]["group-by"], ",")
        keep_columns = get(parsed_args["table"]["aggregate"], "keep-columns", nothing)
        if keep_columns !== nothing
            keep_columns = split(keep_columns, ",")
        end
        count_column = get(parsed_args["table"]["aggregate"], "count-column", "count")
        aggregate_tsv(
            parsed_args["table"]["aggregate"]["input"],
            parsed_args["table"]["aggregate"]["output"];
            group_by=group_by, keep_columns=keep_columns, count_column=count_column
        )
    end

    function handle_unique(parsed_args, _, _)
        columns = split(parsed_args["table"]["unique"]["columns"], ",")
        unique_tsv(parsed_args["table"]["unique"]["input"], parsed_args["table"]["unique"]["output"]; columns=columns)
    end

    function handle_sort(parsed_args, _, _)
        columns = split(parsed_args["table"]["sort"]["columns"], ",")
        reverse = get(parsed_args["table"]["sort"], "reverse", false)
        sort_tsv(parsed_args["table"]["sort"]["input"], parsed_args["table"]["sort"]["output"]; columns=columns, reverse=reverse)
    end

    function handle_filter(parsed_args, _, _)
        column = parsed_args["table"]["filter"]["column"]
        pattern = get(parsed_args["table"]["filter"], "pattern", nothing)
        operator = get(parsed_args["table"]["filter"], "operator", nothing)
        threshold = get(parsed_args["table"]["filter"], "threshold", nothing)
        filter_tsv(parsed_args["table"]["filter"]["input"], parsed_args["table"]["filter"]["output"], column;
                    pattern=pattern, operator=operator, threshold=threshold)
    end

    function handle_select(parsed_args, _, _)
        columns = split(parsed_args["table"]["select"]["columns"], ",")
        select_tsv(parsed_args["table"]["select"]["input"], parsed_args["table"]["select"]["output"]; columns=columns)
    end

    function handle_fasta_export(parsed_args, immunediscover_module, _)
        @info "Exporting TSV to FASTA"
        immunediscover_module.Fasta.extract_sequences_to_fasta(
            parsed_args["table"]["fasta"]["input"],
            parsed_args["table"]["fasta"]["output"];
            colname = parsed_args["table"]["fasta"]["colname"],
            colseq = parsed_args["table"]["fasta"]["colseq"],
            coldesc = parsed_args["table"]["fasta"]["coldesc"],
            filter_pattern = parsed_args["table"]["fasta"]["filter"],
            desc_filter_pattern = parsed_args["table"]["fasta"]["desc-filter"],
            cleanup_pattern = parsed_args["table"]["fasta"]["cleanup"],
            sort_by_name = !parsed_args["table"]["fasta"]["no-sort"],
            mincase = parsed_args["table"]["fasta"]["mincase"],
            case_col = parsed_args["table"]["fasta"]["case-col"],
            unique_sequences = parsed_args["table"]["fasta"]["unique-sequences"]
        )
    end

    function handle_collect(parsed_args, _, _)
        @info "Collecting TSV files"
        pattern = parsed_args["table"]["collect"]["pattern"]
        output = parsed_args["table"]["collect"]["output"]
        files = glob(pattern)
        @info "Found $(length(files)) files matching pattern $pattern"
        if length(files) > 0
            collected = DataFrame[]
            first_file_columns = nothing
            for file in files
                df = CSV.File(file, delim='\t') |> DataFrame
                df[:, :file] .= file
                if first_file_columns === nothing
                    first_file_columns = names(df)
                else
                    @assert first_file_columns == names(df) "Column names in $file do not match first file"
                end
                push!(collected, df)
            end
            collected_df = vcat(collected...)
            CSV.write(output, collected_df, delim='\t', compress=true)
            @info "Collected $(nrow(collected_df)) rows into $output"
        else
            @warn "No files found matching pattern $pattern"
        end
    end

    function handle_exclude(parsed_args, immunediscover_module, always_gz)
        @info "Exclude"
        db = immunediscover_module.load_fasta(parsed_args["table"]["exclude"]["fasta"])
        colname = parsed_args["table"]["exclude"]["colname"]
        colseq = parsed_args["table"]["exclude"]["colseq"]
        data_df = CSV.File(parsed_args["table"]["exclude"]["input"], delim='\t') |> DataFrame
        discard_seqs = Set{String}()
        for row in eachrow(data_df)
            for (name, seq) in db
                if occursin(row[colseq], seq)
                    @info "Allele $(row[colname]) is a substring of $(name)"
                    push!(discard_seqs, row[colseq])
                end
                if occursin(seq, row[colseq])
                    @info "Allele $(name) is a substring of $(row[colname])"
                    push!(discard_seqs, row[colseq])
                end
            end
        end
        for seq in discard_seqs
            @info "Discarding sequences matching: $seq"
            filter!(x -> x[colseq] != seq, data_df)
        end
        output = always_gz(parsed_args["table"]["exclude"]["output"])
        CSV.write(output, data_df, compress=true, delim='\t')
    end

    const TABLE_HANDLERS = Dict{String, Function}(
        "outerjoin"  => handle_outerjoin,
        "leftjoin"   => handle_leftjoin,
        "transform"  => handle_transform,
        "aggregate"  => handle_aggregate,
        "unique"     => handle_unique,
        "sort"       => handle_sort,
        "filter"     => handle_filter,
        "select"     => handle_select,
        "fasta"      => handle_fasta_export,
        "collect"    => handle_collect,
        "exclude"    => handle_exclude,
    )

    function handle_table(parsed_args, immunediscover_module, always_gz)
        subcmd = get(parsed_args["table"], "%COMMAND%", "")
        handler = get(TABLE_HANDLERS, subcmd, nothing)
        if handler !== nothing
            handler(parsed_args, immunediscover_module, always_gz)
        else
            @warn "Unknown table subcommand: $subcmd"
        end
    end
end
