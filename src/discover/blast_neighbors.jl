"Distance between two cores: Hamming when equal length, else Levenshtein."
core_distance(a::AbstractString, b::AbstractString) =
    length(a) == length(b) ? sum(ca != cb for (ca, cb) in zip(a, b); init=0) :
                             Align.core_edit_distance(String(a), String(b))

"""
    neighbor_stats(cores, reads; max_parents) -> (nn_dist, parent_ratio)

For each distinct core, find its nearest *more-abundant* core (its likely "parent"): a
small `nn_dist` with a large `parent_ratio` (= parent reads ÷ this core's reads) marks an
error satellite of a dominant allele. Cores with no more-abundant neighbour get
`nn_dist=-1`, `parent_ratio=1.0`. Only the `max_parents` most-abundant candidates are
examined as potential parents (bounds cost; parents are abundant by definition).
"""
function neighbor_stats(cores::AbstractVector{<:AbstractString}, reads::AbstractVector{<:Integer};
                        max_parents::Int=128)
    n = length(cores)
    nn = fill(-1, n)
    pr = ones(Float64, n)
    order = sortperm(reads, rev=true)            # most abundant first
    for ii in 1:n
        i = order[ii]
        best_d = typemax(Int)
        best_reads = 0
        considered = 0
        for jj in 1:(ii - 1)                      # candidates with ≥ reads
            j = order[jj]
            reads[j] > reads[i] || continue       # parent must be strictly more abundant
            considered += 1
            considered > max_parents && break
            d = core_distance(cores[i], cores[j])
            if d < best_d || (d == best_d && reads[j] > best_reads)
                best_d, best_reads = d, reads[j]
            end
        end
        if best_reads > 0
            nn[i] = best_d
            pr[i] = best_reads / reads[i]
        end
    end
    return nn, pr
end

"""
    satellite_score(nn_dist, parent_ratio) -> Float64

Heuristic satellite score in [0, 1] from proximity to a more-abundant neighbour and read
imbalance. Rises with `parent_ratio` (≥ 20 ≈ full support weight) and nearness (0–1 bp).
"""
function satellite_score(nn_dist::Integer, parent_ratio::Real)
    nn_dist < 0 && return 0.0
    nn_dist > 1 && return 0.0
    pr = Float64(parent_ratio)
    pr < 1.0 && return 0.0
    proximity = 1.0 / (Float64(nn_dist) + 1.0)
    support = clamp(log(pr) / log(20.0), 0.0, 1.0)
    return proximity * support
end

"Derived flag: `satellite_score ≥ 0.5`."
likely_satellite(nn_dist::Integer, parent_ratio::Real) =
    satellite_score(nn_dist, parent_ratio) >= 0.5

"""
    add_neighbor_stats!(df; max_parents)

Add `nn_dist`, `parent_ratio`, `satellite_score`, and `likely_satellite` columns (computed
per gene over the distinct accepted cores; rejected rows keep sentinels -1 / 1.0 / 0 / false).
"""
function add_neighbor_stats!(df::DataFrame; max_parents::Int=128)
    df[!, :nn_dist] = fill(-1, nrow(df))
    df[!, :parent_ratio] = ones(Float64, nrow(df))
    df[!, :satellite_score] = zeros(Float64, nrow(df))
    df[!, :likely_satellite] = falses(nrow(df))
    acc = df[df.reject_reason .== "", :]
    nrow(acc) == 0 && return df
    stats = Dict{String,Tuple{Int,Float64}}()       # core -> (nn_dist, parent_ratio)
    for gdf in groupby(acc, :gene)
        cores = String[]
        reads = Int[]
        seen = Set{String}()
        for r in eachrow(gdf)
            c = String(r.aln_qseq)
            c in seen && continue
            push!(seen, c); push!(cores, c); push!(reads, r.n_reads_total)
        end
        nd, prr = neighbor_stats(cores, reads; max_parents=max_parents)
        for k in eachindex(cores)
            stats[cores[k]] = (nd[k], prr[k])
        end
    end
    @inbounds for i in 1:nrow(df)
        df.reject_reason[i] == "" || continue
        s = get(stats, String(df.aln_qseq[i]), nothing)
        s === nothing && continue
        df.nn_dist[i], df.parent_ratio[i] = s
        df.satellite_score[i] = satellite_score(s[1], s[2])
        df.likely_satellite[i] = likely_satellite(s[1], s[2])
    end
    return df
end
