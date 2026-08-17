# Affix consensus: majority-vote extension from the gene boundary. Prefix/Suffix share
# one AffixSide hierarchy used by both consensus and trimming.

mutable struct AlignmentStats
    total_attempts::Int
    prefix_failures::Int
    suffix_failures::Int
    successful_trims::DefaultDict{String, Vector{String}}
    failed_genes::DefaultDict{String, Vector{String}}
end

AlignmentStats() = AlignmentStats(0, 0, 0,
    DefaultDict{String, Vector{String}}(Vector{String}),
    DefaultDict{String, Vector{String}}(Vector{String}))

"""Merge other into main (for combining thread-local stats after parallel map)."""
function merge_stats!(main::AlignmentStats, other::AlignmentStats)
    main.total_attempts += other.total_attempts
    main.prefix_failures += other.prefix_failures
    main.suffix_failures += other.suffix_failures
    for (k, v) in other.successful_trims
        append!(main.successful_trims[k], v)
    end
    for (k, v) in other.failed_genes
        append!(main.failed_genes[k], v)
    end
    return main
end

abstract type AffixSide end
struct Prefix <: AffixSide end
struct Suffix <: AffixSide end

"""Majority base at one flank column; ties broken lexicographically for stability."""
function majority_base(counts::Dict{Char, Int})
    best_c = first(sort(collect(keys(counts)); by=c -> (-counts[c], c)))
    best_c => counts[best_c]
end

affix_column(a::AbstractString, col::Int, ::Prefix) = a[length(a) - col + 1]
affix_column(a::AbstractString, col::Int, ::Suffix) = a[col]

push_consensus_char!(out::Vector{Char}, c::Char, ::Prefix) = pushfirst!(out, c)
push_consensus_char!(out::Vector{Char}, c::Char, ::Suffix) = push!(out, c)

"""
    consensus_affix(affixes, side; min_fraction) -> String

Shared majority-vote extension from the gene boundary outward.
`Prefix`: 5' flanks, right-aligned, grow leftward.
`Suffix`: 3' flanks, left-aligned, grow rightward.
"""
function consensus_affix(affixes::AbstractVector{<:AbstractString}, side::AffixSide;
                         min_fraction::Real=0.5)
    isempty(affixes) && return ""
    max_len = maximum(length, affixes)
    out = Char[]
    for col in 1:max_len
        counts = Dict{Char, Int}()
        n = 0
        for a in affixes
            length(a) >= col || continue
            c = affix_column(a, col, side)
            counts[c] = get(counts, c, 0) + 1
            n += 1
        end
        n == 0 && break
        best_c, best_n = majority_base(counts)
        best_n / n > min_fraction || break
        push_consensus_char!(out, best_c, side)
    end
    return String(out)
end

consensus_prefix(prefixes::AbstractVector{<:AbstractString}; min_fraction::Real=0.5) =
    consensus_affix(prefixes, Prefix(); min_fraction=min_fraction)

consensus_suffix(suffixes::AbstractVector{<:AbstractString}; min_fraction::Real=0.5) =
    consensus_affix(suffixes, Suffix(); min_fraction=min_fraction)
