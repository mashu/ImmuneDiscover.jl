abstract type CompareOp end
struct Lt <: CompareOp end
struct Le <: CompareOp end
struct Ge <: CompareOp end
struct Gt <: CompareOp end

compare_op(s::AbstractString) = compare_op(Val{Symbol(s)}())
compare_op(::Val{:<}) = Lt()
compare_op(::Val{:<=}) = Le()
compare_op(::Val{:>=}) = Ge()
compare_op(::Val{:>}) = Gt()
compare_op(::Val{S}) where {S} = error("Invalid operator: $S. Must be one of: <, <=, >=, >")

compare(::Lt, x, t) = x < t
compare(::Le, x, t) = x <= t
compare(::Ge, x, t) = x >= t
compare(::Gt, x, t) = x > t

op_label(::Lt) = "<"
op_label(::Le) = "<="
op_label(::Ge) = ">="
op_label(::Gt) = ">"

struct RegexFilter
    pattern::Regex
end
struct NumericFilter{Op<:CompareOp}
    op::Op
    threshold::Float64
end

row_keep(::Missing, ::RegexFilter) = false
row_keep(x, f::RegexFilter) = occursin(f.pattern, string(x))
row_keep(::Missing, ::NumericFilter) = false
row_keep(x::Number, f::NumericFilter) = compare(f.op, x, f.threshold)

as_float(::Missing) = missing
as_float(x::Number) = Float64(x)
as_float(x::AbstractString) = something(tryparse(Float64, x), missing)
as_float(x) = something(tryparse(Float64, string(x)), missing)

unparsed_non_number(::Missing) = false
unparsed_non_number(::Number) = false
unparsed_non_number(_) = true

filter_mask(col, f::RegexFilter) = map(x -> row_keep(x, f), col)
function filter_mask(col, f::NumericFilter)
    numeric = as_float.(col)
    any(ismissing, numeric) && any(unparsed_non_number, col) &&
        error("Column cannot be converted to numeric values for filtering")
    return map(x -> row_keep(x, f), numeric)
end

"""
    filter_tsv(input_file, output_file, column; pattern=absent, operator=absent, threshold=absent)

Filter TSV file by column using regex or numeric operations.
"""
function filter_tsv(input_file, output_file, column; pattern=absent, operator=absent, threshold=absent)
    @info "Loading input file: $input_file"
    df = load_tsv(input_file)
    @info "Input file: $(nrow(df)) rows, $(ncol(df)) columns"
    require_column(df, column)
    spec = filter_spec(optional(pattern), optional(operator), optional(threshold))
    filtered_df = apply_column_filter(df, column, spec)
    @info "Filter results: $(nrow(filtered_df)) rows kept, $(nrow(df) - nrow(filtered_df)) rows removed"
    output_gz = write_tsv(output_file, filtered_df)
    @info "Filtered data saved to: $output_gz"
    return filtered_df
end

filter_spec(::Absent, ::Absent, ::Absent) =
    error("Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering")
filter_spec(p::Present, ::Present, _) =
    error("Cannot use both regex pattern and numeric operator/threshold. Choose one filtering method.")
filter_spec(p::Present, ::Absent, ::Present) =
    error("Cannot use both regex pattern and numeric operator/threshold. Choose one filtering method.")
filter_spec(p::Present, ::Absent, ::Absent) = RegexFilter(Regex(p.value))
filter_spec(::Absent, op::Present, t::Present) = NumericFilter(compare_op(op.value), Float64(t.value))
filter_spec(::Absent, ::Present, ::Absent) =
    error("Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering")
filter_spec(::Absent, ::Absent, ::Present) =
    error("Must provide either --pattern for regex filtering or both --operator and --threshold for numeric filtering")

function apply_column_filter(df, column, f::RegexFilter)
    @info "Filtering column '$column' using regex pattern: $(f.pattern)"
    return df[filter_mask(df[!, column], f), :]
end
function apply_column_filter(df, column, f::NumericFilter)
    @info "Filtering column '$column' using numeric operation: $(op_label(f.op)) $(f.threshold)"
    return df[filter_mask(df[!, column], f), :]
end
