module table
    using CSV
    using DataFrames
    using Logging
    using Glob

    export outerjoin_tsv, leftjoin_tsv, filter_tsv, transform_tsv, aggregate_tsv, unique_tsv, sort_tsv, select_tsv
    export handle_table

    # --- Unified join implementation to eliminate redundancy ---

    """
        join_tsv(left_file, right_file, output_file; how, left_keys, right_keys, left_prefix, right_prefix, left_select, right_select)

    Perform a join (outer or left) on two TSV files by specified key columns.
    Dispatches on `how` to select join type.
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

        @info "Performing $how join"
        result_df = if how === :outer
            outerjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)
        elseif how === :left
            leftjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)
        else
            error("Unsupported join type: $how. Use :outer or :left.")
        end

        @info "Result: $(nrow(result_df)) rows, $(ncol(result_df)) columns"

        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, result_df, compress=true, delim='\t')
        @info "Output saved to: $output_gz"
        return result_df
    end

    # Public API delegates to the unified implementation
    function outerjoin_tsv(left_file, right_file, output_file; kwargs...)
        join_tsv(left_file, right_file, output_file; how=:outer, kwargs...)
    end

    function leftjoin_tsv(left_file, right_file, output_file; kwargs...)
        join_tsv(left_file, right_file, output_file; how=:left, kwargs...)
    end

    # Helper function to get basename without extension
    function splitext(filename)
        if occursin(".", filename)
            parts = split(filename, ".")
            if length(parts) > 1
                return join(parts[1:end-1], "."), parts[end]
            end
        end
        return filename, ""
    end

    """
        is_numeric_column(col)

    Check if a DataFrame column contains numeric data by inspecting its non-missing type.
    Uses dispatch-friendly pattern instead of eltype/isa.
    """
    function is_numeric_column(col)
        for val in col
            ismissing(val) && continue
            return val isa Number
        end
        return false
    end

    """
        parse_numeric_value(val)

    Attempt to parse a value as Float64. Returns nothing on failure (no try-catch).
    """
    function parse_numeric_value(val)
        ismissing(val) && return nothing
        val isa Number && return Float64(val)
        result = tryparse(Float64, string(val))
        return result
    end

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
            # FIX: Use dispatch-friendly numeric detection instead of eltype
            col_data = df[!, column]
            if is_numeric_column(col_data)
                numeric_values = col_data
            else
                # FIX: Use tryparse instead of try-catch for numeric conversion
                parsed = [parse_numeric_value(x) for x in col_data]
                if any(isnothing, parsed)
                    error("Column '$column' cannot be converted to numeric values for filtering")
                end
                numeric_values = [p for p in parsed]
            end

            @info "Filtering column '$column' using numeric operation: $operator $threshold"

            compare_fn = if operator == "<"
                (x, t) -> ismissing(x) ? false : x < t
            elseif operator == "<="
                (x, t) -> ismissing(x) ? false : x <= t
            elseif operator == ">="
                (x, t) -> ismissing(x) ? false : x >= t
            elseif operator == ">"
                (x, t) -> ismissing(x) ? false : x > t
            else
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

    function handle_table(parsed_args, immunediscover_module, always_gz)
        subcmd = get(parsed_args["table"], "%COMMAND%", "")

        if subcmd == "outerjoin"
            @info "Performing outer join"
            _handle_join(parsed_args, "outerjoin", :outer)
        elseif subcmd == "leftjoin"
            @info "Performing left join"
            _handle_join(parsed_args, "leftjoin", :left)
        elseif subcmd == "transform"
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
        elseif subcmd == "aggregate"
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
        elseif subcmd == "unique"
            columns = split(parsed_args["table"]["unique"]["columns"], ",")
            unique_tsv(parsed_args["table"]["unique"]["input"], parsed_args["table"]["unique"]["output"]; columns=columns)
        elseif subcmd == "sort"
            columns = split(parsed_args["table"]["sort"]["columns"], ",")
            reverse = get(parsed_args["table"]["sort"], "reverse", false)
            sort_tsv(parsed_args["table"]["sort"]["input"], parsed_args["table"]["sort"]["output"]; columns=columns, reverse=reverse)
        elseif subcmd == "filter"
            column = parsed_args["table"]["filter"]["column"]
            pattern = get(parsed_args["table"]["filter"], "pattern", nothing)
            operator = get(parsed_args["table"]["filter"], "operator", nothing)
            threshold = get(parsed_args["table"]["filter"], "threshold", nothing)
            filter_tsv(parsed_args["table"]["filter"]["input"], parsed_args["table"]["filter"]["output"], column;
                        pattern=pattern, operator=operator, threshold=threshold)
        elseif subcmd == "select"
            columns = split(parsed_args["table"]["select"]["columns"], ",")
            select_tsv(parsed_args["table"]["select"]["input"], parsed_args["table"]["select"]["output"]; columns=columns)
        elseif subcmd == "fasta"
            @info "Exporting TSV to FASTA"
            immunediscover_module.fasta.extract_sequences_to_fasta(
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
        elseif subcmd == "collect"
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
        elseif subcmd == "exclude"
            @info "Exclude"
            db = immunediscover_module.load_fasta(parsed_args["table"]["exclude"]["fasta"])
            colname = parsed_args["table"]["exclude"]["colname"]
            colseq = parsed_args["table"]["exclude"]["colseq"]
            data_df = CSV.File(parsed_args["table"]["exclude"]["input"], delim='\t') |> DataFrame
            discard = Set{Tuple{String,String,String}}()
            for row in eachrow(data_df)
                for (name, seq) in db
                    if occursin(row[colseq], seq)
                        @info "Allele $(row[colname]) is a substring of $(name)"
                        push!(discard, (name, row[colname], row[colseq]))
                    end
                    if occursin(seq, row[colseq])
                        @info "Allele $(name) is a substring of $(row[colname])"
                        push!(discard, (name, row[colname], row[colseq]))
                    end
                end
            end
            for (name, allele_name, seq) in discard
                @info "Discarding fasta $allele_name with name $name"
                filter!(x -> x[colseq] != seq, data_df)
            end
            output = always_gz(parsed_args["table"]["exclude"]["output"])
            CSV.write(output, data_df, compress=true, delim='\t')
        end
    end

    # Internal helper to avoid duplicating join CLI dispatch logic
    function _handle_join(parsed_args, subcmd, how)
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
end
