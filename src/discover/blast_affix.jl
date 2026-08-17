# Affix consensus (read flanks) and affix trimming (core extraction) share one Prefix/Suffix
# dispatch hierarchy. Previously these were two colliding `AffixSide` definitions in one file.

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

# --- Consensus: grow from the gene boundary outward ---

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

const AFFIX_ANCHOR_K = 12

"""Up to three k-mers per reference (5', middle, 3') for candidate filtering."""
function gene_anchors(ref::AbstractString, k::Int)
    lr = length(ref)
    lr == 0 && return SubString{String}[]
    lr <= k && return [SubString(ref, 1, lr)]
    mid = (lr - k) ÷ 2 + 1
    return unique([SubString(ref, 1, k), SubString(ref, mid, mid + k - 1), SubString(ref, lr - k + 1, lr)])
end

function build_affix_anchor_index(refs::AbstractVector{<:AbstractString}, k::Int)
    index = Dict{SubString{String}, Vector{Int}}()
    for (gi, ref) in enumerate(refs)
        for anchor in gene_anchors(ref, k)
            push!(get!(index, anchor, Int[]), gi)
        end
    end
    return index
end

function affix_candidate_genes(gs::AbstractString, anchor_index::Dict{SubString{String}, Vector{Int}}, k::Int)
    gl = length(gs)
    gl == 0 && return Set{Int}()
    candidates = Set{Int}()
    if gl <= k
        for gi in get(anchor_index, SubString(gs, 1, gl), ())
            push!(candidates, gi)
        end
        return candidates
    end
    @inbounds for pos in 1:(gl - k + 1)
        anchor = SubString(gs, pos, pos + k - 1)
        for gi in get(anchor_index, anchor, ())
            push!(candidates, gi)
        end
    end
    return candidates
end

function affix_hits_for_read(gs::AbstractString, refs::AbstractVector{String},
                             anchor_index::Dict{SubString{String}, Vector{Int}},
                             forward_extension::Int, reverse_extension::Int)
    hits = Tuple{Int, String, String}[]
    gl = length(gs)
    for gi in affix_candidate_genes(gs, anchor_index, AFFIX_ANCHOR_K)
        for m in each_exact_span(refs[gi], gs)
            start_pos, end_pos = extrema(m)
            pre = start_pos > 1 ?
                  String(gs[max(1, start_pos - forward_extension):start_pos - 1]) : ""
            suf = end_pos < gl ?
                  String(gs[end_pos + 1:min(end_pos + reverse_extension, gl)]) : ""
            push!(hits, (gi, pre, suf))
        end
    end
    return hits
end

function accumulate_affixes(db, demux_df; forward_extension=20, reverse_extension=20,
                            min_affix_fraction::Real=0.5)
    names = String[first(p) for p in db]
    refs = String[last(p) for p in db]
    n = length(refs)
    anchor_index = build_affix_anchor_index(refs, AFFIX_ANCHOR_K)
    prefixes = [String[] for _ in 1:n]
    suffixes = [String[] for _ in 1:n]

    prog = Progress(nrow(demux_df); desc="Collecting affixes")
    row_hits = Folds.map(eachrow(demux_df)) do row
        next!(prog)
        affix_hits_for_read(row.genomic_sequence, refs, anchor_index,
                            forward_extension, reverse_extension)
    end
    for hits in row_hits
        for (gi, pre, suf) in hits
            isempty(pre) || push!(prefixes[gi], pre)
            isempty(suf) || push!(suffixes[gi], suf)
        end
    end

    singleton = Vector{Tuple{String, String, String, String}}(undef, n)
    for gi in 1:n
        if isempty(prefixes[gi]) && isempty(suffixes[gi])
            singleton[gi] = (names[gi], refs[gi], "", "")
        else
            common_prefix = consensus_prefix(prefixes[gi]; min_fraction=min_affix_fraction)
            common_suffix = consensus_suffix(suffixes[gi]; min_fraction=min_affix_fraction)
            singleton[gi] = (names[gi], common_prefix * refs[gi] * common_suffix,
                             common_prefix, common_suffix)
        end
    end
    return singleton
end

function save_extended(extended_Ds, fasta_path)
    base_affixes = Vector{Tuple{String, String, String}}()
    open(fasta_path, "w") do io
        for (name, sequence, prefix, suffix) in extended_Ds
            write(io, ">$name\n$sequence\n")
            push!(base_affixes, (name, prefix, suffix))
        end
    end
    return base_affixes
end

# --- Trimming: align affix to the query and keep the opposite side as core ---

side_label(::Prefix) = "prefix"
side_label(::Suffix) = "suffix"

increment_failures!(stats, ::Prefix) = (stats.prefix_failures += 1)
increment_failures!(stats, ::Suffix) = (stats.suffix_failures += 1)

function check_affix_quality_warning(affix_length::Int, quality_threshold::Float64)
    if affix_length > 0 && affix_length <= 20 && quality_threshold > 0.5
        @warn "Quality threshold $(round(quality_threshold * 100, digits=1))% might be too strict for short affixes ($affix_length nt). Consider lowering --minquality."
    end
end

function safe_pairalign(seq1::LongDNA{4}, seq2::LongDNA{4}, scoremodel::AffineGapScoreModel)
    if length(seq1) == 0 || length(seq2) == 0
        @warn "Alignment skipped: one or both sequences are empty"
        return nothing
    end
    return pairalign(SemiGlobalAlignment(), seq1, seq2, scoremodel)
end

function remove_gaps(query::LongDNA{4})::LongDNA{4}
    return LongDNA{4}(filter(nt -> nt != DNA_Gap, query))
end

"""For prefix: core is everything after the last aligned affix position."""
function extract_core_range(aligned_affix, aligned_query, ::Prefix)
    boundary = findlast(x -> x != DNA_Gap, aligned_affix)
    boundary === nothing && return nothing, "No match found"
    start = boundary + 1
    start > length(aligned_query) && return nothing, "Query too short after trimming"
    return aligned_query[start:end], ""
end

"""For suffix: core is everything before the first aligned affix position."""
function extract_core_range(aligned_affix, aligned_query, ::Suffix)
    boundary = findfirst(x -> x != DNA_Gap, aligned_affix)
    boundary === nothing && return nothing, "No match found"
    boundary <= 1 && return nothing, "Starts too early"
    return aligned_query[1:boundary-1], ""
end

"""
    trim_one_affix(affix, query, stats, scoremodel, side; min_quality, sseqid)

Align `affix` to `query`, check quality, and extract the core sequence on the
opposite side. Returns the trimmed query or nothing on failure.
"""
function trim_one_affix(affix::LongDNA{4}, query::LongDNA{4}, stats, scoremodel, side::AffixSide;
                        min_quality=0.75, sseqid="")
    length(affix) == 0 && return query

    label = side_label(side)
    aln_result = safe_pairalign(affix, query, scoremodel)
    if aln_result === nothing
        push!(stats.failed_genes[sseqid], "$(titlecase(label)) alignment failed")
        increment_failures!(stats, side)
        return nothing
    end

    pairs = collect(alignment(aln_result))
    positions = findall(p -> first(p) != DNA_Gap, pairs)
    if isempty(positions)
        push!(stats.failed_genes[sseqid], "No $label content in alignment")
        increment_failures!(stats, side)
        return nothing
    end

    matches = sum(first(pairs[i]) == last(pairs[i]) for i in positions)
    quality = matches / length(positions)
    if quality < min_quality
        push!(stats.failed_genes[sseqid], "Poor $label alignment quality ($(round(quality * 100, digits=1))% match)")
        increment_failures!(stats, side)
        return nothing
    end

    aligned_affix = first.(pairs)
    aligned_query = last.(pairs)
    core, msg = extract_core_range(aligned_affix, aligned_query, side)
    if core === nothing
        push!(stats.failed_genes[sseqid], "$(titlecase(label)): $msg")
        increment_failures!(stats, side)
        return nothing
    end

    result = remove_gaps(LongDNA{4}(join(core)))
    if length(result) == 0
        push!(stats.failed_genes[sseqid], "Empty core after $label trimming")
        increment_failures!(stats, side)
        return nothing
    end
    return result
end

function trim_sequence(query::LongDNA{4}, prefix::LongDNA{4}, suffix::LongDNA{4}, stats,
    scoremodel::AffineGapScoreModel=AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1);
    min_quality=0.75, sseqid="")

    if length(query) == 0
        push!(stats.failed_genes[sseqid], "Empty query sequence")
        return nothing
    end
    stats.total_attempts += 1

    partial = trim_one_affix(prefix, query, stats, scoremodel, Prefix();
                             min_quality=min_quality, sseqid=sseqid)
    partial === nothing && return nothing

    partial = trim_one_affix(suffix, partial, stats, scoremodel, Suffix();
                             min_quality=min_quality, sseqid=sseqid)
    partial === nothing && return nothing

    if length(partial) < 10
        push!(stats.failed_genes[sseqid], "Core sequence too short ($(length(partial)) < 10 nt)")
        stats.prefix_failures += 1
        return nothing
    end

    push!(stats.successful_trims[sseqid], "Trimmed successfully")
    return partial
end

function trim_and_align_sequence(query::String, prefix::String, suffix::String, reference::String, stats; min_quality=0.75, sseqid="")
    query_dna = BioSequences.LongDNA{4}(query)
    prefix_dna = BioSequences.LongDNA{4}(prefix)
    suffix_dna = BioSequences.LongDNA{4}(suffix)
    trimmed_dna = trim_sequence(query_dna, prefix_dna, suffix_dna, stats, min_quality=min_quality, sseqid=sseqid)
    trimmed_dna === nothing && return "", -1
    trimmed = replace(String(trimmed_dna), '-' => "")
    distance = Align.core_edit_distance(trimmed, reference)
    return trimmed, distance
end
