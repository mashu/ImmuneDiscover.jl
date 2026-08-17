# ========================== Per-mode match schemas ==========================
# The output schema follows the gene and mode the user is searching, so a V search never
# carries D columns. Each mode/gene produces a CONCRETE row type, which keeps the per-read
# hot loop type-stable (a `Vector{<concrete NamedTuple>}`, not `Vector{NamedTuple}`):
#   - extension mode (any gene): ExtRow — prefix/sequence/suffix + their lengths;
#   - RSS mode: the gene-specific flank NamedTuple from `extract_flanking` (V: prefix +
#     heptamer/spacer/nonamer; J: suffix + heptamer/spacer/nonamer; D: pre_*/post_*).
# The `--rss` selection is applied as a column PROJECTION (`project_rss!`) rather than by
# varying the per-row type, so grouping/output keep exactly the requested RSS elements.

const ExtRow = NamedTuple{(:well, :case, :db_name, :prefix, :sequence, :suffix, :prefix_len, :suffix_len),
                          Tuple{String, String, String, String, String, String, Int, Int}}

"Keep only the gene-appropriate RSS columns for `--rss` (V/J: mandatory flank + selected RSS; D: both sides)."
function project_rss!(df::DataFrame, ::VGene, rss)
    optional = intersect(["heptamer", "spacer", "nonamer"], rss)
    select!(df, vcat(["well", "case", "db_name", "prefix", "sequence"], optional))
end
function project_rss!(df::DataFrame, ::JGene, rss)
    optional = intersect(["heptamer", "spacer", "nonamer"], rss)
    select!(df, vcat(["well", "case", "db_name", "sequence", "suffix"], optional))
end
function d_rss_columns(rss)
    cols = String[]
    "heptamer" in rss && append!(cols, ["pre_heptamer", "post_heptamer"])
    "spacer" in rss && append!(cols, ["pre_spacer", "post_spacer"])
    "nonamer" in rss && append!(cols, ["pre_nonamer", "post_nonamer"])
    for c in ("pre_heptamer", "pre_spacer", "pre_nonamer", "post_heptamer", "post_spacer", "post_nonamer")
        c in rss && !(c in cols) && push!(cols, c)
    end
    return cols
end
function project_rss!(df::DataFrame, ::DGene, rss)
    select!(df, vcat(["well", "case", "db_name", "sequence"], d_rss_columns(rss)))
end


# ========================== Merge helpers ==========================

function merge_counts!(dst::Dict{Tuple{String,String},Int}, src::Dict{Tuple{String,String},Int})
    for (k, v) in src
        dst[k] = get(dst, k, 0) + v
    end
    return dst
end

# ========================== Range validation ==========================

function validate_range(start_pos, end_pos, seq_length)
    (start_pos < 1 || end_pos > seq_length || start_pos > end_pos) && error("Invalid range")
end

# ========================== extract_flanking — RSS mode (dispatch) ==========================

extract_flanking(gs::String, range::Tuple{Int,Int}, gt::GeneType, n::Int) =
    extract_flanking(gs, range, gt, n, nothing)

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::VGene, n::Int, ::Nothing)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
    return (prefix=prefix, sequence=seq,
            heptamer=flank_slice(gs, ep+1, ep+7),
            spacer=flank_slice(gs, ep+8, ep+30),
            nonamer=flank_slice(gs, ep+31, ep+39))
end

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::JGene, n::Int, ::Nothing)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    return (nonamer=flank_slice(gs, sp-39, sp-31),
            spacer=flank_slice(gs, sp-30, sp-8),
            heptamer=flank_slice(gs, sp-7, sp-1),
            suffix=flank_slice(gs, ep+1, ep+n),
            sequence=seq)
end

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ::Nothing)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    return (pre_nonamer=flank_slice(gs, sp-28, sp-20),
            pre_spacer=flank_slice(gs, sp-19, sp-8),
            pre_heptamer=flank_slice(gs, sp-7, sp-1),
            sequence=seq,
            post_heptamer=flank_slice(gs, ep+1, ep+7),
            post_spacer=flank_slice(gs, ep+8, ep+19),
            post_nonamer=flank_slice(gs, ep+20, ep+28))
end

# ========================== extract_flanking — Extension mode (dispatch) ==========================

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::VGene, n::Int, ext::Int)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
    return (prefix=prefix, sequence=seq, suffix=gs[ep+1:min(ep+ext,length(gs))])
end

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::JGene, n::Int, ext::Int)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    return (prefix=flank_slice(gs, sp-ext, sp-1), sequence=seq, suffix=flank_slice(gs, ep+1, ep+n))
end

function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ext::Int)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    return (prefix=flank_slice(gs, sp-ext, sp-1), sequence=seq, suffix=flank_slice(gs, ep+1, ep+ext))
end

# ========================== extract_flanking — Per-side extension (dispatch) ==========================

function extract_flanking(gs::String, range::Tuple{Int,Int}, gt::GeneType, n::Int, ext::Int, prefE::Int, sufE::Int)
    sp, ep = range; validate_range(sp, ep, length(gs))
    seq = gs[sp:ep]
    return extract_perside(gs, sp, ep, seq, gt, n, prefE, sufE)
end

function extract_perside(gs, sp, ep, seq, ::VGene, n, prefE, sufE)
    prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
    return (prefix=prefix, sequence=seq, suffix=gs[ep+1:min(ep+sufE,length(gs))])
end
function extract_perside(gs, sp, ep, seq, ::JGene, n, prefE, sufE)
    return (prefix=flank_slice(gs, sp-prefE, sp-1), sequence=seq, suffix=flank_slice(gs, ep+1, ep+n))
end
function extract_perside(gs, sp, ep, seq, ::DGene, n, prefE, sufE)
    return (prefix=flank_slice(gs, sp-prefE, sp-1), sequence=seq, suffix=flank_slice(gs, ep+1, ep+sufE))
end

# ========================== extension_overlaps_border (dispatch) ==========================

function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::VGene, ext::Int, border::Int)
    (border <= 0 || ext <= 0) && return false
    rbs = max(1, rl - border + 1)
    ee = min(ep + ext, rl)
    return (ep + 1 <= ee) && (ee >= rbs)
end

function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::JGene, ext::Int, border::Int)
    (border <= 0 || ext <= 0) && return false
    lbe = min(border, rl)
    es = max(1, sp - ext)
    return (es <= sp - 1) && (es <= lbe)
end

function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::DGene, ext::Int, border::Int)
    (border <= 0 || ext <= 0) && return false
    lbe = min(border, rl); rbs = max(1, rl - border + 1)
    les = max(1, sp - ext); ree = min(ep + ext, rl)
    return ((les <= sp-1) && (les <= lbe)) || ((ep+1 <= ree) && (ree >= rbs))
end

# ========================== Calibration helpers (dispatch) ==========================

function collect_allowed_lengths!(pre, suf, ::VGene, gb, ext, sp, ep, rl, border)
    rbs = max(1, rl - border + 1)
    push!(get!(suf, gb, Int[]), min(ext, max(0, rbs - 1 - ep)))
end
function collect_allowed_lengths!(pre, suf, ::JGene, gb, ext, sp, ep, rl, border)
    lbe = min(border, rl)
    push!(get!(pre, gb, Int[]), min(ext, max(0, sp - 1 - lbe)))
end
function collect_allowed_lengths!(pre, suf, ::DGene, gb, ext, sp, ep, rl, border)
    lbe = min(border, rl); rbs = max(1, rl - border + 1)
    push!(get!(pre, gb, Int[]), min(ext, max(0, sp - 1 - lbe)))
    push!(get!(suf, gb, Int[]), min(ext, max(0, rbs - 1 - ep)))
end

function safe_quantile(values::AbstractVector{<:Integer}, p::Float64, fallback::Int)
    isempty(values) && return fallback
    sorted = sort(values)
    idx = max(1, min(length(sorted), ceil(Int, (1 - p) * length(sorted))))
    return sorted[idx]
end

function assign_calibrated!(pgp, pgs, ::VGene, pre_v, suf_v, ext, ap)
    for (g, arr) in suf_v; pgs[g] = safe_quantile(arr, ap, ext); end
end
function assign_calibrated!(pgp, pgs, ::JGene, pre_v, suf_v, ext, ap)
    for (g, arr) in pre_v; pgp[g] = safe_quantile(arr, ap, ext); end
end
function assign_calibrated!(pgp, pgs, ::DGene, pre_v, suf_v, ext, ap)
    for (g, arr) in pre_v; pgp[g] = safe_quantile(arr, ap, ext); end
    for (g, arr) in suf_v; pgs[g] = safe_quantile(arr, ap, ext); end
end

# ========================== Border rejection (dispatch) ==========================

function should_reject_border(::VGene, pgp, pgs, gb, ext, sp, ep, rl, border)
    rbs = max(1, rl - border + 1)
    return min(ep + get(pgs, gb, ext), rl) >= rbs
end
function should_reject_border(::JGene, pgp, pgs, gb, ext, sp, ep, rl, border)
    lbe = min(border, rl)
    return max(1, sp - get(pgp, gb, ext)) <= lbe
end
function should_reject_border(::DGene, pgp, pgs, gb, ext, sp, ep, rl, border)
    lbe = min(border, rl); rbs = max(1, rl - border + 1)
    return (max(1, sp - get(pgp, gb, ext)) <= lbe) || (min(ep + get(pgs, gb, ext), rl) >= rbs)
end

# ========================== Ratio utilities ==========================

# Ratio with an explicit 0-denominator convention: control genes (0 reference count/median)
# yield Inf so they pass the downstream min-ratio filters by design, never NaN.
safe_ratio(num, den) = den == 0 ? Inf : num / den
