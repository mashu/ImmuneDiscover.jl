module Filters

using DataFrames

export FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative, CustomFilter
export GermlineFilter, passes, apply_filters!, add_group_ratio!
export init_rejection_columns!, mark_rejected!, annotate_rejections!, accepted

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
    MinThreshold(:full_count, 5, "Min full cluster count (--minfullcount)"),
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

# --- Annotate path: record WHY a row would be dropped instead of dropping it ---------------
# Used by the discovery pipelines to emit a full table (every candidate + the reason it was
# rejected) alongside the filtered table (rows with an empty reason). The first rejection a
# row hits wins, so reasons accumulate across stages without overwriting an earlier one.

"""
    init_rejection_columns!(df; reason_col=:reject_reason, stage_col=:reject_stage)

Ensure the string `reason_col`/`stage_col` columns exist (default ""). Returns `df`.
"""
function init_rejection_columns!(df::DataFrame; reason_col::Symbol=:reject_reason, stage_col::Symbol=:reject_stage)
    reason_col in propertynames(df) || (df[!, reason_col] = fill("", nrow(df)))
    stage_col in propertynames(df) || (df[!, stage_col] = fill("", nrow(df)))
    return df
end

"""
    mark_rejected!(df, fail_mask, reason, stage; reason_col=:reject_reason, stage_col=:reject_stage)

Mark rows where `fail_mask` is true AND not already rejected with `reason`/`stage`. Returns `df`.
"""
function mark_rejected!(df::DataFrame, fail_mask::AbstractVector{Bool}, reason::AbstractString, stage::AbstractString;
                        reason_col::Symbol=:reject_reason, stage_col::Symbol=:reject_stage)
    init_rejection_columns!(df; reason_col=reason_col, stage_col=stage_col)
    reasons = df[!, reason_col]; stages = df[!, stage_col]
    @inbounds for i in eachindex(fail_mask)
        if fail_mask[i] && isempty(reasons[i])
            reasons[i] = reason
            stages[i] = stage
        end
    end
    return df
end

"""
    annotate_rejections!(df, criteria; stage="filter", reason_col=:reject_reason, stage_col=:reject_stage)

Mark each not-yet-rejected row with the label of the FIRST `criteria` it fails (or leave it
"" if it passes all). Like `GermlineFilter` but it annotates rather than removes. Returns `df`.
"""
function annotate_rejections!(df::DataFrame, criteria::Vector{<:FilterCriterion};
                              stage::AbstractString="filter",
                              reason_col::Symbol=:reject_reason, stage_col::Symbol=:reject_stage)
    init_rejection_columns!(df; reason_col=reason_col, stage_col=stage_col)
    reasons = df[!, reason_col]; stages = df[!, stage_col]
    for (i, row) in enumerate(eachrow(df))
        isempty(reasons[i]) || continue
        for criterion in criteria
            if !passes(row, criterion)
                reasons[i] = criterion.label
                stages[i] = stage
                break
            end
        end
    end
    return df
end

"""
    accepted(df; reason_col=:reject_reason) -> view of rows that passed every stage.
"""
accepted(df::DataFrame; reason_col::Symbol=:reject_reason) =
    filter(row -> isempty(getproperty(row, reason_col)), df)

end
