function aggregate_tsv(input_file, output_file; group_by, keep_columns=absent, count_column="count")
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
    require_columns(df, group_by)
    keep = keep_column_list(optional(keep_columns), df, group_by)
    require_columns(df, keep)
    selected_columns = vcat(group_by, keep)
    df_selected = df[:, selected_columns]
    agg_specs = Any[nrow => count_column]
    for col in keep
        (col in group_by || col == count_column) && continue
        push!(agg_specs, col => first => col)
    end
    grouped_df = combine(groupby(df_selected, group_by), agg_specs...)
    @info "Aggregated from $(nrow(df)) rows to $(nrow(grouped_df)) unique groups"
    write_tsv(output_file, grouped_df)
    return grouped_df
end

keep_column_list(::Absent, df, group_by) = [col for col in names(df) if !(col in group_by)]
keep_column_list(s::Present{<:AbstractVector}, _, _) = String.(s.value)
keep_column_list(s::Present{<:AbstractString}, _, _) = parse_csv_list(s.value)

function unique_tsv(input_file, output_file; columns)
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    require_columns(df, columns)
    unique_df = unique(df[:, columns])
    @info "Unique rows: $(nrow(unique_df)) (removed $(nrow(df) - nrow(unique_df)) duplicates)"
    write_tsv(output_file, unique_df)
    return unique_df
end

function sort_tsv(input_file, output_file; columns, reverse=false)
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    require_columns(df, columns)
    sort!(df, columns, rev=reverse)
    @info "Sorted $(nrow(df)) rows by $columns ($(reverse ? "desc" : "asc"))"
    write_tsv(output_file, df)
    return df
end

function select_tsv(input_file, output_file; columns)
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    require_columns(df, columns)
    selected_df = df[:, columns]
    @info "Selected $(ncol(selected_df)) columns from $(ncol(df))"
    write_tsv(output_file, selected_df)
    return selected_df
end
