# Collect 5'/3' read flanks around exact gene hits and vote a consensus affix per reference.

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

read_affix(gs, span, extension, ::Prefix) =
    flank_slice(gs, first(span) - extension, first(span) - 1)
read_affix(gs, span, extension, ::Suffix) =
    flank_slice(gs, last(span) + 1, last(span) + extension)

function affix_hits_for_read(gs::AbstractString, refs::AbstractVector{String},
                             anchor_index::Dict{SubString{String}, Vector{Int}},
                             forward_extension::Int, reverse_extension::Int)
    hits = Tuple{Int, String, String}[]
    for gi in affix_candidate_genes(gs, anchor_index, AFFIX_ANCHOR_K)
        for span in each_exact_span(refs[gi], gs)
            pre = read_affix(gs, span, forward_extension, Prefix())
            suf = read_affix(gs, span, reverse_extension, Suffix())
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
        singleton[gi] = extended_record(names[gi], refs[gi], prefixes[gi], suffixes[gi];
                                        min_affix_fraction=min_affix_fraction)
    end
    return singleton
end

function extended_record(name, ref, prefixes, suffixes; min_affix_fraction)
    isempty(prefixes) && isempty(suffixes) && return (name, ref, "", "")
    common_prefix = consensus_prefix(prefixes; min_fraction=min_affix_fraction)
    common_suffix = consensus_suffix(suffixes; min_fraction=min_affix_fraction)
    return (name, common_prefix * ref * common_suffix, common_prefix, common_suffix)
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
