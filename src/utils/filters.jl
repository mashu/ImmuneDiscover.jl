module Filters

using DataFrames
using Logging

export FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative, CustomFilter
export GermlineFilter, passes, apply_filters!

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
    for criterion in gf.criteria
        before = nrow(df)
        filter!(row -> passes(row, criterion), df)
        @info "$(criterion.label): $(nrow(df))/$before rows kept"
    end
    return df
end

"""
    apply_filters!(df, criteria...)

Convenience: apply individual FilterCriterion values without constructing a GermlineFilter.
"""
function apply_filters!(df::DataFrame, criteria::FilterCriterion...)
    GermlineFilter(collect(FilterCriterion, criteria))(df)
end

end
