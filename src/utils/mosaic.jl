module Mosaic
    using DataFrames
    using ..Align: core_edit_distance

    export chimera_score, refs_by_gene, add_chimera_scores!, gene_from_allele

    const MIN_CHIMERA_SEGMENT = 4

    "Gene prefix from an allele name (`IGHV1-2*01` → `IGHV1-2`)."
    gene_from_allele(name::AbstractString) = first(split(strip(String(name)), '*'))

    "Edit cost of `seg` anchored to the start (`:start`) or end (`:end`) of `ref`."
    function segment_anchor_cost(seg::AbstractString, ref::AbstractString, side::Symbol)
        Ls, Lr = length(seg), length(ref)
        Ls == 0 && return 0
        L = min(Ls, Lr)
        if side === :start
            tail = Ls - L
            return core_edit_distance(seg[1:L], ref[1:L]) + tail
        end
        side === :end || error("side must be :start or :end")
        tail = Ls - L
        return core_edit_distance(seg[Ls - L + 1:Ls], ref[Lr - L + 1:Lr]) + tail
    end

    """
        chimera_score(seq, refs) -> Float64

    Heuristic mosaic score in [0, 1]: excess benefit of explaining `seq` as a left prefix of
    one reference plus a right suffix of another, vs a single reference. `refs` is a vector of
    sequences or `(name, seq)` tuples (same gene).
    """
    function chimera_score(seq::AbstractString, refs::AbstractVector{<:AbstractString})
        L = length(seq)
        L < 2 * MIN_CHIMERA_SEGMENT && return 0.0
        length(refs) < 2 && return 0.0
        single = minimum(core_edit_distance(seq, ref) for ref in refs)
        single == 0 && return 0.0
        best_mosaic = typemax(Int)
        for split in MIN_CHIMERA_SEGMENT:(L - MIN_CHIMERA_SEGMENT)
            left = seq[1:split]
            right = seq[split+1:L]
            left_cost = minimum(segment_anchor_cost(left, ref, :start) for ref in refs)
            right_cost = minimum(segment_anchor_cost(right, ref, :end) for ref in refs)
            best_mosaic = min(best_mosaic, left_cost + right_cost)
        end
        excess = single - best_mosaic
        excess <= 0 && return 0.0
        return clamp(excess / single, 0.0, 1.0)
    end

    chimera_score(seq::AbstractString, refs::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}}) =
        chimera_score(seq, [String(s) for (_, s) in refs])

    "Group FASTA `(name, seq)` pairs by gene for per-gene chimera scoring."
    function refs_by_gene(db_seqs)
        d = Dict{String, Vector{Tuple{String, String}}}()
        for (n, s) in db_seqs
            g = gene_from_allele(n)
            push!(get!(d, g, Tuple{String, String}[]), (strip(String(n)), String(s)))
        end
        return d
    end

    """
        add_chimera_scores!(df, db_by_gene; seq_col, gene_col)

    Append `chimera_score` per row, caching by `(gene, sequence)`.
    """
    function add_chimera_scores!(df::DataFrame, db_by_gene::AbstractDict;
                                 seq_col::Symbol, gene_col::Symbol)
        df[!, :chimera_score] = zeros(Float64, nrow(df))
        cache = Dict{Tuple{String, String}, Float64}()
        @inbounds for i in 1:nrow(df)
            gene = String(df[i, gene_col])
            seq = String(df[i, seq_col])
            key = (gene, seq)
            score = get(cache, key, nothing)
            if score === nothing
                refs = get(db_by_gene, gene, nothing)
                score = refs === nothing || length(refs) < 2 ? 0.0 : chimera_score(seq, refs)
                cache[key] = score
            end
            df.chimera_score[i] = score
        end
        return df
    end
end
