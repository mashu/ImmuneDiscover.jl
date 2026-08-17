function transform_tsv(input_file, output_file; column, pattern, replacement, new_column=absent)
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
    require_column(df, column)
    @info "Transforming column '$column' using pattern: $pattern"
    regex = Regex(pattern)
    original_values = df[!, column]
    transformed_values = String[]
    extracted_values = String[]
    dest = optional(new_column)
    for value in original_values
        m = match(regex, string(value))
        apply_match!(transformed_values, extracted_values, m, replacement, dest, string(value))
    end
    changed_count = sum(string.(original_values) .!= transformed_values)
    @info "Transformed $changed_count values in column '$column'"
    df[!, column] = transformed_values
    add_extracted_column!(df, dest, extracted_values)
    output_gz = write_tsv(output_file, df)
    @info "Transformed data saved to: $output_gz"
    return df
end

apply_match!(transformed, extracted, ::Nothing, _, ::Absent, original) =
    (push!(transformed, original); push!(extracted, ""))
apply_match!(transformed, extracted, ::Nothing, _, ::Present, original) =
    (push!(transformed, original); push!(extracted, ""))
function apply_match!(transformed, extracted, m::RegexMatch, replacement, dest, _)
    result = replacement
    for (i, capture) in enumerate(m.captures)
        capture === nothing && continue
        result = replace(result, "\\$i" => capture)
    end
    push!(transformed, result)
    push!(extracted, extracted_text(dest, m))
end

extracted_text(::Absent, _) = ""
extracted_text(::Present, m::RegexMatch) = join(c for c in m.captures if c !== nothing)

add_extracted_column!(df, ::Absent, _) = df
function add_extracted_column!(df, col::Present, extracted)
    df[!, col.value] = extracted
    @info "Added new column '$(col.value)'"
    return df
end
