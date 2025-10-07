module table
    using CSV
    using DataFrames
    using Logging

    export outerjoin_tsv, leftjoin_tsv, filter_tsv, transform_tsv, aggregate_tsv, unique_tsv, sort_tsv, select_tsv

    """
        outerjoin_tsv(left_file, right_file, output_file; left_keys, right_keys, left_prefix, right_prefix, left_select, right_select)

    Perform outer join on two TSV files by specified key columns.
    
    # Arguments
    - `left_file`: Path to the left TSV file
    - `right_file`: Path to the right TSV file  
    - `output_file`: Path for the output TSV file
    - `left_keys`: Vector of column names to join on from left file
    - `right_keys`: Vector of column names to join on from right file (defaults to left_keys)
    - `left_prefix`: Optional prefix for left file columns (defaults to basename of left_file)
    - `right_prefix`: Optional prefix for right file columns (defaults to basename of right_file)
    - `left_select`: Optional vector of column names to select from left file (defaults to all columns)
    - `right_select`: Optional vector of column names to select from right file (defaults to all columns)
    """
    function outerjoin_tsv(left_file, right_file, output_file; 
                           left_keys, right_keys=nothing, left_prefix=nothing, right_prefix=nothing, 
                           left_select=nothing, right_select=nothing)
        
        @info "Loading left file: $left_file"
        left_df = CSV.File(left_file, delim='\t') |> DataFrame
        
        @info "Loading right file: $right_file"
        right_df = CSV.File(right_file, delim='\t') |> DataFrame
        
        # Use left_keys as right_keys if not specified
        if right_keys === nothing
            right_keys = left_keys
        end
        
        # Set default prefixes to nothing (no prefixing by default)
        if left_prefix === nothing
            left_prefix = ""
        end
        if right_prefix === nothing
            right_prefix = ""
        end
        
        # Select columns if specified
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
        @info "Prefixes - Left: $left_prefix, Right: $right_prefix"
        
        # Add prefixes to column names (except join keys) only if prefix is not empty
        left_df_prefixed = copy(left_df)
        right_df_prefixed = copy(right_df)
        
        for col in names(left_df)
            if !(col in left_keys) && !isempty(left_prefix)
                new_name = "$(left_prefix)_$(col)"
                rename!(left_df_prefixed, col => new_name)
            end
        end
        
        for col in names(right_df)
            if !(col in right_keys) && !isempty(right_prefix)
                new_name = "$(right_prefix)_$(col)"
                rename!(right_df_prefixed, col => new_name)
            end
        end
        
        # Perform outer join
        @info "Performing outer join"
        result_df = outerjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)
        
        @info "Result: $(nrow(result_df)) rows, $(ncol(result_df)) columns"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, result_df, compress=true, delim='\t')
        @info "Output saved to: $output_gz"
        
        return result_df
    end

    """
        leftjoin_tsv(left_file, right_file, output_file; left_keys, right_keys, left_prefix, right_prefix, left_select, right_select)

    Perform left join on two TSV files by specified key columns.
    
    # Arguments
    - `left_file`: Path to the left TSV file
    - `right_file`: Path to the right TSV file  
    - `output_file`: Path for the output TSV file
    - `left_keys`: Vector of column names to join on from left file
    - `right_keys`: Vector of column names to join on from right file (defaults to left_keys)
    - `left_prefix`: Optional prefix for left file columns (defaults to basename of left_file)
    - `right_prefix`: Optional prefix for right file columns (defaults to basename of right_file)
    - `left_select`: Optional vector of column names to select from left file (defaults to all columns)
    - `right_select`: Optional vector of column names to select from right file (defaults to all columns)
    """
    function leftjoin_tsv(left_file, right_file, output_file; 
                          left_keys, right_keys=nothing, left_prefix=nothing, right_prefix=nothing, 
                          left_select=nothing, right_select=nothing)
        
        @info "Loading left file: $left_file"
        left_df = CSV.File(left_file, delim='\t') |> DataFrame
        
        @info "Loading right file: $right_file"
        right_df = CSV.File(right_file, delim='\t') |> DataFrame
        
        # Use left_keys as right_keys if not specified
        if right_keys === nothing
            right_keys = left_keys
        end
        
        # Set default prefixes to nothing (no prefixing by default)
        if left_prefix === nothing
            left_prefix = ""
        end
        if right_prefix === nothing
            right_prefix = ""
        end
        
        # Select columns if specified
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
        @info "Prefixes - Left: $left_prefix, Right: $right_prefix"
        
        # Add prefixes to column names (except join keys) only if prefix is not empty
        left_df_prefixed = copy(left_df)
        right_df_prefixed = copy(right_df)
        
        for col in names(left_df)
            if !(col in left_keys) && !isempty(left_prefix)
                new_name = "$(left_prefix)_$(col)"
                rename!(left_df_prefixed, col => new_name)
            end
        end
        
        for col in names(right_df)
            if !(col in right_keys) && !isempty(right_prefix)
                new_name = "$(right_prefix)_$(col)"
                rename!(right_df_prefixed, col => new_name)
            end
        end
        
        # Perform left join
        @info "Performing left join"
        result_df = leftjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)
        
        @info "Result: $(nrow(result_df)) rows, $(ncol(result_df)) columns"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, result_df, compress=true, delim='\t')
        @info "Output saved to: $output_gz"
        
        return result_df
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
        filter_tsv(input_file, output_file, column; pattern=nothing, operator=nothing, threshold=nothing)

    Filter TSV file by column using regex or numeric operations.
    Mode is automatically detected based on provided parameters:
    - If `pattern` is provided: regex mode (string filtering)
    - If `operator` and `threshold` are provided: numeric mode (numeric filtering)
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `column`: Name of the column to filter on
    - `pattern`: Regex pattern for string filtering (optional)
    - `operator`: Numeric operator - "<", "<=", ">=", ">" (optional)
    - `threshold`: Numeric threshold value (optional)
    """
    function filter_tsv(input_file, output_file, column; pattern=nothing, operator=nothing, threshold=nothing)
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if column exists
        if !(column in names(df))
            error("Column '$column' not found in input file. Available columns: $(join(names(df), ", "))")
        end
        
        # Auto-detect mode based on provided parameters
        if pattern !== nothing && (operator !== nothing || threshold !== nothing)
            error("Cannot use both regex pattern and numeric operator/threshold. Choose one filtering method.")
        elseif pattern !== nothing
            # Regex mode
            @info "Filtering column '$column' using regex pattern: $pattern"
            # Compile regex once before the loop for performance
            regex = Regex(pattern)
            # Handle missing values by skipping them (treating them as not matching the pattern)
            filter_mask = map(x -> ismissing(x) ? false : occursin(regex, x), df[!, column])
        elseif operator !== nothing && threshold !== nothing
            # Numeric mode
            # Check if column can be converted to numeric
            numeric_values = try
                # Try to convert to numeric, handling both string and already numeric columns
                if eltype(df[!, column]) <: Number
                    df[!, column]
                else
                    parse.(Float64, string.(df[!, column]))
                end
            catch
                error("Column '$column' cannot be converted to numeric values for filtering")
            end
            
            @info "Filtering column '$column' using numeric operation: $operator $threshold"
            
            # Apply numeric filter based on operator, handling missing values
            if operator == "<"
                filter_mask = map(x -> ismissing(x) ? false : x < threshold, numeric_values)
            elseif operator == "<="
                filter_mask = map(x -> ismissing(x) ? false : x <= threshold, numeric_values)
            elseif operator == ">="
                filter_mask = map(x -> ismissing(x) ? false : x >= threshold, numeric_values)
            elseif operator == ">"
                filter_mask = map(x -> ismissing(x) ? false : x > threshold, numeric_values)
            else
                error("Invalid operator: $operator. Must be one of: <, <=, >=, >")
            end
        else
            error("Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering")
        end
        
        # Apply filter
        filtered_df = df[filter_mask, :]
        
        original_count = nrow(df)
        filtered_count = nrow(filtered_df)
        removed_count = original_count - filtered_count
        
        @info "Filter results: $filtered_count rows kept, $removed_count rows removed"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, filtered_df, compress=true, delim='\t')
        @info "Filtered data saved to: $output_gz"
        
        return filtered_df
    end

    """
        transform_tsv(input_file, output_file, column; pattern, replacement, new_column=nothing)

    Transform a TSV file by applying regex pattern with capture groups to replace matched strings.
    Optionally add a new column with captured groups.
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `column`: Name of the column to transform
    - `pattern`: Regex pattern with capture groups (e.g., "ID_(\\d+)_(\\w+)" to capture ID_123_ABC)
    - `replacement`: Replacement string using capture groups (e.g., "\\1-\\2" to get "123-ABC")
    - `new_column`: Optional name for new column containing captured groups (e.g., "extracted_id")
    """
    function transform_tsv(input_file, output_file; column, pattern, replacement, new_column=nothing)
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if column exists
        if !(column in names(df))
            error("Column '$column' not found in input file. Available columns: $(join(names(df), ", "))")
        end
        
        @info "Transforming column '$column' using pattern: $pattern"
        @info "Replacement pattern: $replacement"
        
        # Apply regex transformation
        original_values = df[!, column]
        @info "Julia replacement pattern: $replacement"
        
        # Manual regex replacement with capture groups
        transformed_values = String[]
        extracted_values = String[]
        regex = Regex(pattern)
        
        for value in original_values
            m = match(regex, string(value))
            if m !== nothing
                # Replace \1, \2, etc. with actual capture groups
                result = replacement
                for i in 1:length(m.captures)
                    if m.captures[i] !== nothing
                        result = replace(result, "\\$i" => m.captures[i])
                    end
                end
                push!(transformed_values, result)
                
                # Extract captured groups for new column if requested
                if new_column !== nothing
                    extracted = join([m.captures[i] for i in 1:length(m.captures) if m.captures[i] !== nothing], "")
                    push!(extracted_values, extracted)
                else
                    push!(extracted_values, "")
                end
            else
                push!(transformed_values, string(value))
                push!(extracted_values, "")
            end
        end
        
        # Count how many values were actually changed
        changed_count = sum(original_values .!= transformed_values)
        @info "Transformed $changed_count values in column '$column'"
        
        # Update the dataframe
        df[!, column] = transformed_values
        
        # Add new column with extracted values if requested
        if new_column !== nothing
            df[!, new_column] = extracted_values
            @info "Added new column '$new_column' with extracted values"
        end
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, df, compress=true, delim='\t')
        @info "Transformed data saved to: $output_gz"
        
        return df
    end

    """
        aggregate_tsv(input_file, output_file; group_by, keep_columns=nothing, count_column="count")

    Aggregate unique rows by selected columns, optionally keeping specific columns and adding a count.
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `group_by`: Vector of column names to group by
    - `keep_columns`: Optional vector of additional columns to keep (defaults to all non-group columns)
    - `count_column`: Name for the count column (defaults to "count")
    """
    function aggregate_tsv(input_file, output_file; group_by, keep_columns=nothing, count_column="count")
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if group_by columns exist
        for col in group_by
            if !(col in names(df))
                error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
            end
        end
        
        # Determine which columns to keep
        if keep_columns === nothing
            # Keep all columns that are not in group_by
            keep_columns = [col for col in names(df) if !(col in group_by)]
        else
            # Check if keep_columns exist
            for col in keep_columns
                if !(col in names(df))
                    error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
                end
            end
        end
        
        @info "Grouping by: $group_by"
        @info "Keeping columns: $keep_columns"
        
        # Select only the columns we need
        selected_columns = vcat(group_by, keep_columns)
        df_selected = df[:, selected_columns]
        
        # Group by the specified columns and count
        grouped_df = combine(groupby(df_selected, group_by), nrow => count_column)
        
        # Add any additional columns that were kept (taking first value from each group)
        if !isempty(keep_columns)
            for col in keep_columns
                if !(col in names(grouped_df))
                    first_values = combine(groupby(df_selected, group_by), col => first => col)
                    grouped_df = leftjoin(grouped_df, first_values, on=group_by)
                end
            end
        end
        
        original_count = nrow(df)
        aggregated_count = nrow(grouped_df)
        @info "Aggregated from $original_count rows to $aggregated_count unique groups"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, grouped_df, compress=true, delim='\t')
        @info "Aggregated data saved to: $output_gz"
        
        return grouped_df
    end

    """
        unique_tsv(input_file, output_file; columns)

    Create unique rows by selecting specific columns (like SELECT DISTINCT).
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `columns`: Vector of column names to select for uniqueness
    """
    function unique_tsv(input_file, output_file; columns)
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if columns exist
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
            end
        end
        
        @info "Selecting unique rows based on columns: $columns"
        
        # Select only the specified columns and get unique rows
        unique_df = unique(df[:, columns])
        
        original_count = nrow(df)
        unique_count = nrow(unique_df)
        removed_count = original_count - unique_count
        
        @info "Unique rows: $unique_count rows (removed $removed_count duplicates)"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, unique_df, compress=true, delim='\t')
        @info "Unique data saved to: $output_gz"
        
        return unique_df
    end

    """
        sort_tsv(input_file, output_file; columns, reverse=false)

    Sort a TSV file by specified columns.
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `columns`: Vector of column names to sort by (in order of priority)
    - `reverse`: Whether to sort in descending order (defaults to false for ascending)
    """
    function sort_tsv(input_file, output_file; columns, reverse=false)
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if columns exist
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
            end
        end
        
        @info "Sorting by columns: $columns"
        @info "Sort order: $(reverse ? "descending" : "ascending")"
        
        # Sort the dataframe
        sort!(df, columns, rev=reverse)
        
        @info "Sorted $(nrow(df)) rows"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, df, compress=true, delim='\t')
        @info "Sorted data saved to: $output_gz"
        
        return df
    end

    """
        select_tsv(input_file, output_file; columns)

    Select specific columns from a TSV file and save to output file.
    
    # Arguments
    - `input_file`: Path to the input TSV file
    - `output_file`: Path for the output TSV file
    - `columns`: Vector of column names to select
    """
    function select_tsv(input_file, output_file; columns)
        
        @info "Loading input file: $input_file"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
        
        # Check if all specified columns exist
        for col in columns
            if !(col in names(df))
                error("Column '$col' not found in input file. Available columns: $(join(names(df), ", "))")
            end
        end
        
        @info "Selecting columns: $columns"
        
        # Select only the specified columns
        selected_df = df[:, columns]
        
        @info "Selected $(ncol(selected_df)) columns from $(ncol(df)) original columns"
        
        # Save result
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, selected_df, compress=true, delim='\t')
        @info "Selected data saved to: $output_gz"
        
        return selected_df
    end
end
