module table
    using CSV
    using DataFrames
    using Logging

    export outerjoin_tsv, leftjoin_tsv, filter_tsv

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
        result_df = outerjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys)
        
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
        result_df = leftjoin(left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys)
        
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
            filter_mask = occursin.(Regex(pattern), df[!, column])
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
            
            # Apply numeric filter based on operator
            if operator == "<"
                filter_mask = numeric_values .< threshold
            elseif operator == "<="
                filter_mask = numeric_values .<= threshold
            elseif operator == ">="
                filter_mask = numeric_values .>= threshold
            elseif operator == ">"
                filter_mask = numeric_values .> threshold
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
end
