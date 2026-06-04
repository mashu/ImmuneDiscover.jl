module Filters

using DataFrames

export FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative, CustomFilter
export GermlineFilter, passes, apply_filters!, add_group_ratio!

abstract type FilterCriterion end

struct MinThreshold <: FilterCriterion
    column::Symbol
    value::Float64
    label::String
end

struct MaxThreshold <: FilterCriterion
    column::Symbol
    value::Float64
    label::String
end

struct MinStringLength <: FilterCriterion
    column::Symbol
    value::Int
    label::String
end

struct NonNegative <: FilterCriterion
    column::Symbol
    label::String
end

struct CustomFilter <: FilterCriterion
    predicate::Function
    label::String
end

passes(row, f::MinThreshold) = getproperty(row, f.column) >= f.value
passes(row, f::MaxThreshold) = getproperty(row, f.column) <= f.value
passes(row, f::MinStringLength) = length(getproperty(row, f.column)) >= f.value
passes(row, f::NonNegative) = getproperty(row, f.column) >= 0
passes(row, f::CustomFilter) = f.predicate(row)

"""
    GermlineFilter(criteria)

Composable filter for germline discovery pipelines. Applies a sequence of
dispatch-based criteria to a DataFrame, logging kept/total counts at each step.

# Usage
```julia
gf = GermlineFilter([
    MinThreshold(:full_count, 5, "Min cluster size"),
    MaxThreshold(:mismatch, 10, "Max edit distance"),
    MinStringLength(:qseq, 290, "Min read length"),
])
gf(df)  # filters df in-place, returns df
```
"""
struct GermlineFilter
    criteria::Vector{FilterCriterion}
end

function (gf::GermlineFilter)(df::DataFrame)
    isempty(gf.criteria) && return df
    start_rows = nrow(df)
    for criterion in gf.criteria
        before = nrow(df)
        filter!(row -> passes(row, criterion), df)
        report_filter_step(criterion.label, nrow(df), before)
    end
    report_filter_summary(start_rows, nrow(df))
    return df
end

# --- Colored diagnostics for the filter cascade (respects the terminal's color support) ---

function report_filter_step(label::AbstractString, kept::Int, before::Int)
    removed = before - kept
    printstyled("  ✓ "; color = :green, bold = true)
    print(rpad(label, 34), " ")
    printstyled(string(kept); color = :cyan, bold = true)
    print("/$before kept")
    removed > 0 && printstyled("  −$removed"; color = :light_red)
    println()
end

function report_filter_summary(start_rows::Int, kept::Int)
    printstyled("  Σ "; color = :blue, bold = true)
    print("kept ")
    printstyled(string(kept); color = (kept == 0 ? :light_red : :green), bold = true)
    print("/$start_rows")
    start_rows > 0 && print(" (", round(100 * kept / start_rows; digits = 1), "%)")
    println()
end

# Compact and pretty `show` for the GermlineFilter functor.
function Base.show(io::IO, gf::GermlineFilter)
    n = length(gf.criteria)
    printstyled(io, "GermlineFilter"; color = :cyan, bold = true)
    print(io, "(", n, n == 1 ? " criterion)" : " criteria)")
end

function Base.show(io::IO, ::MIME"text/plain", gf::GermlineFilter)
    printstyled(io, "GermlineFilter"; color = :cyan, bold = true)
    println(io, " with $(length(gf.criteria)) criteria:")
    for c in gf.criteria
        printstyled(io, "  • "; color = :green)
        printstyled(io, c.label; color = :yellow)
        println(io)
    end
end

"""
    apply_filters!(df, criteria...)

Convenience: apply individual FilterCriterion values without constructing a GermlineFilter.
"""
function apply_filters!(df::DataFrame, criteria::FilterCriterion...)
    GermlineFilter(collect(FilterCriterion, criteria))(df)
end

"""
    add_group_ratio!(df, value_col, group_cols, ratio_col)

Add `ratio_col` = `value_col` divided by its per-group maximum (groups defined by
`group_cols`). This is the standard "allelic ratio within gene" computed across the
pipeline before applying a ratio threshold; factored here so every caller is consistent.
Returns `df`.
"""
function add_group_ratio!(df::DataFrame, value_col::Symbol, group_cols, ratio_col::Symbol)
    transform!(groupby(df, group_cols), value_col => (x -> x ./ maximum(x)) => ratio_col)
    return df
end

end
