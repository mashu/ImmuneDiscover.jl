# ========================== Extension calibration (stage) ==========================

"""
    calibrate_extension(table, query, gt, extension, border, adjust_percent) -> (prefix, suffix)

Per-gene safe prefix/suffix extension lengths so that ≥ `adjust_percent` of reads avoid the
read border. Only used in extension mode with `border > 0` and per-gene adjustment enabled.
"""
function calibrate_extension(table, query, gt::GeneType, extension::Int, border::Int, adjust_percent::Float64)
    per_gene_prefix = Dict{String,Int}()
    per_gene_suffix = Dict{String,Int}()
    @info "Calibrating per-gene extension targeting ≥ $(Int(round(adjust_percent*100)))% safe reads"
    p_cal = Progress(nrow(table))
    tmp = Folds.map(eachrow(table)) do row
        next!(p_cal)
        local pre = Dict{String,Vector{Int}}()
        local suf = Dict{String,Vector{Int}}()
        rl = length(row.genomic_sequence)
        @inbounds for (name, seq) in query
            for m in each_exact_span(seq, row.genomic_sequence)
                gb = first(split(string(name), '*'))
                collect_allowed_lengths!(pre, suf, gt, gb, extension, minimum(m), maximum(m), rl, border)
            end
        end
        (pre=pre, suf=suf)
    end
    pre_vals = Dict{String,Vector{Int}}()
    suf_vals = Dict{String,Vector{Int}}()
    for r in tmp
        for (g, v) in r.pre; haskey(pre_vals, g) ? append!(pre_vals[g], v) : (pre_vals[g] = copy(v)); end
        for (g, v) in r.suf; haskey(suf_vals, g) ? append!(suf_vals[g], v) : (suf_vals[g] = copy(v)); end
    end
    assign_calibrated!(per_gene_prefix, per_gene_suffix, gt, pre_vals, suf_vals, extension, adjust_percent)
    @info "Per-gene extension calibrated for $(length(union(collect(keys(per_gene_prefix)), collect(keys(per_gene_suffix))))) genes"
    return per_gene_prefix, per_gene_suffix
end

# ========================== Matching (stage) ==========================

"""
    collect_matches(table, query, gt, affix, rss, extension, border, adjust, pgp, pgs)
        -> (result_df, totals_all, accepted_all)

Scan every read for exact occurrences of each query allele and emit one match row per accepted
occurrence, in the gene/mode-appropriate schema. Dispatches to a type-stable per-mode kernel
(RSS vs extension); `totals_all`/`accepted_all` hold border-filter tallies (extension only).
"""
function collect_matches(table, query, gt::GeneType, affix::Int, rss, ::Absent,
                         border::Int, adjust::Bool, per_gene_prefix, per_gene_suffix)
    result_df = collect_rss(table, query, gt, affix)
    isempty(result_df) || project_rss!(result_df, gt, rss)
    empty = Dict{Tuple{String,String},Int}()
    return result_df, empty, copy(empty)
end

function collect_matches(table, query, gt::GeneType, affix::Int, rss, e::Present,
                         border::Int, adjust::Bool, per_gene_prefix, per_gene_suffix)
    return collect_extension(table, query, gt, affix, e.value, border, adjust,
                             per_gene_prefix, per_gene_suffix)
end

"RSS mode (no extension): full gene-specific flank rows; `gt` barrier keeps the row type concrete."
function collect_rss(table, query, gt::G, affix::Int) where {G<:GeneType}
    p = Progress(nrow(table))
    per_read = Folds.map(eachrow(table)) do row
        next!(p)
        read_matches_rss(row, query, gt, affix)
    end
    valid = [v for v in per_read if !isempty(v)]
    return isempty(valid) ? DataFrame() : DataFrame(reduce(vcat, valid))
end

"Match rows for one read in RSS mode. The comprehension infers a concrete element type under `G`."
function read_matches_rss(row, query, gt::G, affix::Int) where {G<:GeneType}
    well = string(row.well); case = string(row.case); gs = row.genomic_sequence
    return [merge((well=well, case=case, db_name=string(name)),
                  extract_flanking(gs, (minimum(m), maximum(m)), gt, affix, absent))
            for (name, seq) in query for m in each_exact_span(seq, gs)]
end

"Extension mode: fixed `ExtRow` schema (prefix/sequence/suffix + lengths) plus border tallies."
function collect_extension(table, query, gt::G, affix::Int, extension::Int, border::Int,
                           adjust::Bool, per_gene_prefix, per_gene_suffix) where {G<:GeneType}
    p = Progress(nrow(table))
    result = Folds.map(eachrow(table)) do row
        next!(p)
        well = string(row.well); case = string(row.case); gs = row.genomic_sequence; rl = length(gs)
        matches = ExtRow[]
        totals = Dict{Tuple{String,String},Int}()
        accepted_counts = Dict{Tuple{String,String},Int}()
        @inbounds for (name, seq) in query
            for m in each_exact_span(seq, gs)
                sp = minimum(m); ep = maximum(m)
                gb = first(split(string(name), '*'))
                if border > 0
                    key = (case, gb)
                    totals[key] = get(totals, key, 0) + 1
                    rejected = adjust ?
                        should_reject_border(gt, per_gene_prefix, per_gene_suffix, gb, extension, sp, ep, rl, border) :
                        extension_overlaps_border(sp, ep, rl, gt, extension, border)
                    rejected && continue
                    accepted_counts[key] = get(accepted_counts, key, 0) + 1
                end
                flanks = (adjust && border > 0) ?
                    extract_flanking(gs, (sp, ep), gt, affix, extension, get(per_gene_prefix, gb, extension), get(per_gene_suffix, gb, extension)) :
                    extract_flanking(gs, (sp, ep), gt, affix, extension)
                push!(matches, (well=well, case=case, db_name=string(name),
                                prefix=String(flanks.prefix), sequence=String(flanks.sequence),
                                suffix=String(flanks.suffix),
                                prefix_len=length(flanks.prefix), suffix_len=length(flanks.suffix)))
            end
        end
        (matches=matches, totals=totals, accepted=accepted_counts)
    end

    valid_match_lists = [r.matches for r in result if !isempty(r.matches)]
    result_df = isempty(valid_match_lists) ? DataFrame() : DataFrame(reduce(vcat, valid_match_lists))

    totals_all = Dict{Tuple{String,String},Int}()
    accepted_all = Dict{Tuple{String,String},Int}()
    for r in result
        !isempty(r.totals) && merge_counts!(totals_all, r.totals)
        !isempty(r.accepted) && merge_counts!(accepted_all, r.accepted)
    end
    return result_df, totals_all, accepted_all
end

# ========================== Border statistics (stage) ==========================

"""
    summarize_border_stats!(totals_all, accepted_all, active)

Log the border-filter rejection rate. No process-global state.
"""
function summarize_border_stats!(totals_all, accepted_all, active::Bool)
    (!active || isempty(totals_all)) && return
    tm = 0; tr = 0
    for (c, g) in union(collect(keys(totals_all)), collect(keys(accepted_all)))
        mt = get(totals_all, (c, g), 0)
        a = get(accepted_all, (c, g), 0)
        tm += mt
        tr += mt - a
    end
    @info "Border filter rejected $tr of $tm potential matches ($(tm>0 ? round(100*tr/tm;digits=2) : 0.0)%)"
    return
end
