module Exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds
    using FASTX
    using Statistics
    using ..Filters: FilterCriterion, MinThreshold, MinStringLength, CustomFilter, add_group_ratio!,
                     init_rejection_columns!, mark_rejected!, accepted, passes
    using ..Mosaic: refs_by_gene, add_chimera_scores!
    using ..Data: barplot_if_available, boxplot_if_available, round_floats!
    using ..Report: section, stage_report, report_rejections,
                    filter_quality_report, rss_consistency

    # ========================== GeneType dispatch hierarchy ==========================

    abstract type GeneType end
    struct VGene <: GeneType end
    struct DGene <: GeneType end
    struct JGene <: GeneType end

    const GENE_TYPE_MAP = Dict{String, GeneType}("V" => VGene(), "J" => JGene(), "D" => DGene())

    """
        parse_gene_type(s) -> GeneType

    Convert a gene type string ("V", "J", "D") to a dispatch-ready type.
    """
    function parse_gene_type(s::AbstractString)
        gt = get(GENE_TYPE_MAP, String(s), nothing)
        gt === nothing && error("Invalid gene type: $s")
        return gt
    end

    gene_string(::VGene) = "V"
    gene_string(::JGene) = "J"
    gene_string(::DGene) = "D"

    const GENE_CHAR_MAP = (('V', VGene()), ('D', DGene()), ('J', JGene()))

    """
        gene_type_from_name(name) -> GeneType or nothing

    Infer gene type from an allele/gene name (e.g. "IGHV1-2" → VGene()).
    Returns nothing if no V/D/J is found.
    """
    function gene_type_from_name(name::AbstractString)
        for (ch, gt) in GENE_CHAR_MAP
            occursin(ch, name) && return gt
        end
        return nothing
    end

    export GeneType, VGene, DGene, JGene, parse_gene_type, gene_type_from_name

    # ========================== Typed rows for border statistics ==========================

    const BORDER_ROW = NamedTuple{(:case, :gene, :matched_total, :accepted_total, :rejected_border, :rejected_ratio), Tuple{String, String, Int, Int, Int, Float64}}
    const BORDER_GENE_ROW = NamedTuple{(:gene, :mean_rejected_ratio, :matched_total, :rejected_border, :num_donors), Tuple{String, Float64, Int, Int, Int}}

    const LAST_BORDER_STATS = Ref(Vector{BORDER_ROW}())
    const LAST_BORDER_GENE_STATS = Ref(Vector{BORDER_GENE_ROW}())

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

    "Keep only the gene-appropriate RSS columns for `--rss` (V/J: mandatory flank + selected RSS; D: all)."
    function project_rss!(df::DataFrame, ::VGene, rss)
        optional = intersect(["heptamer", "spacer", "nonamer"], rss)
        select!(df, vcat(["well", "case", "db_name", "prefix", "sequence"], optional))
    end
    function project_rss!(df::DataFrame, ::JGene, rss)
        optional = intersect(["heptamer", "spacer", "nonamer"], rss)
        select!(df, vcat(["well", "case", "db_name", "sequence", "suffix"], optional))
    end
    project_rss!(df::DataFrame, ::DGene, rss) = df   # D keeps its full pre_/post_ RSS context


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

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::VGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
        return (prefix=prefix, sequence=seq,
                heptamer=gs[ep+1:min(ep+7,length(gs))],
                spacer=gs[ep+8:min(ep+30,length(gs))],
                nonamer=gs[ep+31:min(ep+39,length(gs))])
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::JGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (nonamer=gs[max(1,sp-39):max(1,sp-31)],
                spacer=gs[max(1,sp-30):max(1,sp-8)],
                heptamer=gs[max(1,sp-7):(sp-1)],
                suffix=gs[ep+1:min(ep+n,length(gs))],
                sequence=seq)
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (pre_nonamer=gs[max(1,sp-28):max(1,sp-20)],
                pre_spacer=gs[max(1,sp-19):max(1,sp-8)],
                pre_heptamer=gs[max(1,sp-7):(sp-1)],
                sequence=seq,
                post_heptamer=gs[ep+1:min(ep+7,length(gs))],
                post_spacer=gs[ep+8:min(ep+19,length(gs))],
                post_nonamer=gs[ep+20:min(ep+28,length(gs))])
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
        return (prefix=gs[max(1,sp-ext):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+n,length(gs))])
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ext::Int)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (prefix=gs[max(1,sp-ext):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+ext,length(gs))])
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
        return (prefix=gs[max(1,sp-prefE):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+n,length(gs))])
    end
    function extract_perside(gs, sp, ep, seq, ::DGene, n, prefE, sufE)
        return (prefix=gs[max(1,sp-prefE):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+sufE,length(gs))])
    end

    # ========================== String interface (used by tests and CLI) ==========================

    function extract_flanking(gs::String, range::Tuple{Int,Int}, gene_type::String, n::Int, extension::Union{Int,Nothing}=nothing)
        extract_flanking(gs, range, parse_gene_type(gene_type), n, extension)
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, gene_type::String, n::Int,
                              extension::Union{Int,Nothing}, prefix_ext::Union{Int,Nothing}, suffix_ext::Union{Int,Nothing})
        gt = parse_gene_type(gene_type)
        extension === nothing && return extract_flanking(gs, range, gt, n, nothing)
        left_ext = prefix_ext === nothing ? extension : prefix_ext
        right_ext = suffix_ext === nothing ? extension : suffix_ext
        return extract_flanking(gs, range, gt, n, extension, left_ext, right_ext)
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

    function extension_overlaps_border(sp::Int, ep::Int, rl::Int, gene_type::String, ext::Int, border::Int)
        extension_overlaps_border(sp, ep, rl, parse_gene_type(gene_type), ext, border)
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

    function get_ratio(expect_dict, row, ratio)
        row.db_name in keys(expect_dict) && (@info "Skipping allelic ratio filters for $(row.db_name) in case $(row.case)"; return 0.0)
        row.gene in keys(expect_dict) && (@info "Skipping gene ratio filters for $(row.db_name) in case $(row.case)"; return 0.0)
        return ratio
    end

    # Ratio with an explicit 0-denominator convention: control genes (0 reference count/median)
    # yield Inf so they pass the downstream min-ratio filters by design, never NaN.
    safe_ratio(num, den) = den == 0 ? Inf : num / den

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
                m = findfirst(seq, row.genomic_sequence)
                if m !== nothing
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
    function collect_matches(table, query, gt::GeneType, affix::Int, rss, extension,
                             border::Int, adjust::Bool, per_gene_prefix, per_gene_suffix)
        if extension === nothing
            result_df = collect_rss(table, query, gt, affix)
            isempty(result_df) || project_rss!(result_df, gt, rss)
            empty = Dict{Tuple{String,String},Int}()
            return result_df, empty, copy(empty)
        end
        return collect_extension(table, query, gt, affix, extension, border, adjust,
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
                      extract_flanking(gs, (minimum(m), maximum(m)), gt, affix, nothing))
                for (name, seq) in query for m in (findfirst(seq, gs),) if m !== nothing]
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
                m = findfirst(seq, gs)
                m === nothing && continue
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

    Populate `LAST_BORDER_STATS` / `LAST_BORDER_GENE_STATS` (per case×gene and per gene) from the
    border-filter tallies, and log the overall rejection rate. A no-op (empties the refs) when the
    border filter was not active.
    """
    function summarize_border_stats!(totals_all, accepted_all, active::Bool)
        if !active || isempty(totals_all)
            LAST_BORDER_STATS[] = BORDER_ROW[]; LAST_BORDER_GENE_STATS[] = BORDER_GENE_ROW[]
            return
        end
        stats_rows = BORDER_ROW[]
        tm = 0; tr = 0
        for (c, g) in union(collect(keys(totals_all)), collect(keys(accepted_all)))
            mt = get(totals_all, (c, g), 0); a = get(accepted_all, (c, g), 0); r = mt - a
            push!(stats_rows, (case=c, gene=g, matched_total=mt, accepted_total=a, rejected_border=r, rejected_ratio=(mt > 0 ? r/mt : 0.0)))
            tm += mt; tr += r
        end
        @info "Border filter rejected $tr of $tm potential matches ($(tm>0 ? round(100*tr/tm;digits=2) : 0.0)%)"
        LAST_BORDER_STATS[] = stats_rows
        sr = Dict{String,Float64}(); sm = Dict{String,Int}(); srej = Dict{String,Int}(); nd = Dict{String,Int}()
        for row in stats_rows
            g = row.gene; sr[g] = get(sr, g, 0.0) + row.rejected_ratio; sm[g] = get(sm, g, 0) + row.matched_total
            srej[g] = get(srej, g, 0) + row.rejected_border; nd[g] = get(nd, g, 0) + 1
        end
        LAST_BORDER_GENE_STATS[] = [
            (gene=g, mean_rejected_ratio=(nd[g]>0 ? sr[g]/nd[g] : 0.0), matched_total=sm[g], rejected_border=srej[g], num_donors=nd[g])
            for g in keys(sr)]
        return
    end

    # ========================== Counting & metrics (stages) ==========================

    """
        add_counts!(result_df, sequence_lookup) -> df

    Add `full_count` (identical full rows), `count` (per well/case/db_name/sequence), `gene`,
    optional `isin_db`, and the within-gene allelic ratios `full_ratio` / `ratio`.
    """
    function add_counts!(result_df::DataFrame, sequence_lookup)
        df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
        transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
        transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)
        if sequence_lookup !== nothing
            @info "Adding isin_db column based on reference FASTA"
            df[!, :isin_db] = map(row -> get(sequence_lookup, row.sequence, false) ? "" : "Novel", eachrow(df))
        end
        add_group_ratio!(df, :full_count, [:well, :case, :gene], :full_ratio)
        add_group_ratio!(df, :count, [:well, :case, :gene], :ratio)
        return df
    end

    """
        add_quality_metrics!(udf) -> udf

    Add per-candidate-core metrics over distinct rows: `n_donors` (cross-donor recurrence),
    `n_reads_total` (read support), `max_full_ratio` (peak per-donor allelic ratio). Must be
    called on the de-duplicated table so `full_count` is summed once per distinct row.
    """
    function add_quality_metrics!(udf::DataFrame)
        transform!(groupby(udf, :sequence), :case => (x -> length(unique(x))) => :n_donors)
        transform!(groupby(udf, :sequence), :full_count => sum => :n_reads_total)
        transform!(groupby(udf, :sequence), :full_ratio => maximum => :max_full_ratio)
        return udf
    end

    """
        locus_group_stat!(df, groupcols, srccol, destcol, locus, statfn; default=0)

    Within each group, set `destcol` to `statfn` over `srccol` for the rows that both start with
    `locus` and are still accepted (`reject_reason == ""`); `default` when none qualify. Control
    genes (outside `locus`) and already-rejected rows are excluded from the statistic.
    """
    function locus_group_stat!(df::DataFrame, groupcols, srccol::Symbol, destcol::Symbol, locus, statfn; default=0)
        transform!(groupby(df, groupcols)) do g
            fg = filter(r -> startswith(r.db_name, locus) && isempty(r.reject_reason), g)
            DataFrame(destcol => fill(isempty(fg) ? default : statfn(fg[!, srccol]), nrow(g)))
        end
        return df
    end

    """
        add_frequency_columns!(df, locus) -> df

    Add the locus-frequency columns used by the reference-frequency filters: `gene_count` /
    `case_count` (per-well/case totals), `allele_cohort_median` / `gene_cohort_median` (across-donor
    medians), `gene_case_freq` (gene-usage fraction), `allele_cohort_fold` / `gene_cohort_fold`
    (fold-change vs cohort median), and `allelic_ratio` (within-gene allelic ratio). All aggregates
    are over accepted (count/ratio-passing) rows. Requires `reject_reason` (count/ratio annotated first).
    """
    function add_frequency_columns!(df::DataFrame, locus::AbstractString)
        # Per-well/case totals (over accepted, in-locus rows).
        locus_group_stat!(df, [:well, :case, :gene], :count, :gene_count, locus, sum)
        locus_group_stat!(df, [:well, :case], :count, :case_count, locus, sum)
        # Cohort (across-donor) medians — robust central tendency used for the fold-change below.
        locus_group_stat!(df, [:db_name], :count, :allele_cohort_median, locus, median; default=0.0)
        locus_group_stat!(df, [:gene], :gene_count, :gene_cohort_median, locus, median; default=0.0)

        # gene_case_freq = gene-usage fraction in the case (low ⇒ possible deletion; see --deletion).
        df[:, :gene_case_freq] = safe_ratio.(df.gene_count, df.case_count)
        # *_cohort_fold = this donor's count ÷ the allele's/gene's cohort-median (robust fold-change
        # vs typical); a tiny fold flags a sporadic low-support observation.
        df[:, :allele_cohort_fold] = safe_ratio.(df.count, df.allele_cohort_median)
        df[:, :gene_cohort_fold] = safe_ratio.(df.gene_count, df.gene_cohort_median)
        # allelic_ratio = a row's reads as a fraction of its gene's accepted reads in that
        # well+case (the standard within-gene allele-calling signal).
        transform!(groupby(df, [:well, :case, :gene])) do g
            denom = sum((r.count for r in eachrow(g) if isempty(r.reject_reason)); init=0)
            DataFrame(allelic_ratio = safe_ratio.(g.count, denom))
        end
        return df
    end

    # ========================== Filter assembly & annotation (stages) ==========================

    """
        exact_filter_criteria(; mincount, minratio, expect_dict, min_recurrence, min_seqlen, min_peak_ratio)

    The count/ratio criteria (plus optional quality floors) applied to exact candidates. `count`
    is intentionally omitted: it is ≥ `full_count` by construction, so `full_count ≥ mincount`
    already subsumes it. The ratio criterion requires both the full-row and collapsed allelic
    ratios to clear the (control-gene-relaxed) threshold.
    """
    function exact_filter_criteria(; mincount, minratio, expect_dict,
                                   min_recurrence::Int=0, min_seqlen::Int=0, min_peak_ratio::Float64=0.0)
        crit = FilterCriterion[
            MinThreshold(:full_count, Float64(mincount), "min count (--mincount $mincount)"),
            CustomFilter(r -> (thr = get_ratio(expect_dict, r, minratio); r.full_ratio >= thr && r.ratio >= thr),
                         "min allelic ratio (--minratio $minratio)"),
        ]
        min_recurrence > 0 && push!(crit,
            MinThreshold(:n_donors, Float64(min_recurrence), "min donor recurrence (--min-recurrence $min_recurrence)"))
        min_seqlen > 0 && push!(crit,
            MinStringLength(:sequence, min_seqlen, "min sequence length (--min-seqlen $min_seqlen)"))
        min_peak_ratio > 0 && push!(crit,
            MinThreshold(:max_full_ratio, min_peak_ratio, "min peak allelic ratio (--min-peak-ratio $min_peak_ratio)"))
        return crit
    end

    """
        exact_frequency_criteria(mod, expect_dict, deletion_dict, min_allele_fold, min_gene_fold)

    The reference-frequency criteria: allele/gene-case frequency floors (control-gene-aware via
    `get_ratio_threshold`) and the cross-case median ratios.
    """
    function exact_frequency_criteria(mod, expect_dict, deletion_dict, min_allele_fold, min_gene_fold)
        return FilterCriterion[
            CustomFilter(x -> x.allelic_ratio >= mod.get_ratio_threshold(expect_dict, x, type="allelic_ratio"), "allelic ratio"),
            CustomFilter(x -> x.gene_case_freq >= mod.get_ratio_threshold(deletion_dict, x, type="gene_case_freq"), "gene-case frequency"),
            MinThreshold(:allele_cohort_fold, min_allele_fold, "min allele cohort fold (--min-allele-cohort-fold)"),
            MinThreshold(:gene_cohort_fold, min_gene_fold, "min gene cohort fold (--min-gene-cohort-fold)"),
        ]
    end

    """
        annotate_stage!(df, criteria, stage)

    Apply each criterion in turn, marking (not dropping) the first rejection per row and printing
    a per-criterion kept/removed line. Shared shape with the discovery pipelines.
    """
    function annotate_stage!(df::DataFrame, criteria, stage::AbstractString)
        for criterion in criteria
            before = count(isempty, df.reject_reason)
            fail = Bool[!passes(row, criterion) for row in eachrow(df)]
            mark_rejected!(df, fail, criterion.label, stage)
            stage_report(criterion.label, count(isempty, df.reject_reason), before)
        end
        return df
    end

    # ========================== exact_search (orchestrator) ==========================

    """
        exact_search(table, query, gene; kwargs...) -> DataFrame

    Find exact occurrences of each `query` allele in the reads and return the UNFILTERED candidate
    table (one row per distinct flank/sequence, capped at `N` flank records per allele) with counts,
    allelic ratios and quality metrics. Callers annotate/filter as needed (`handle_exact` runs the
    full transparency cascade; `hsmm` keeps `full_count ≥ mincount`).
    """
    function exact_search(table, query, gene; affix=13, rss=["heptamer", "spacer", "nonamer"],
                          extension=nothing, N=10, raw=nothing, sequence_lookup=nothing,
                          border::Int=0, adjust_per_gene_extension::Bool=false, adjust_percent::Float64=1.0)
        gt = parse_gene_type(gene)
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"

        per_gene_prefix = Dict{String,Int}()
        per_gene_suffix = Dict{String,Int}()
        if extension !== nothing && border > 0 && adjust_per_gene_extension
            per_gene_prefix, per_gene_suffix = calibrate_extension(table, query, gt, extension, border, adjust_percent)
        end

        result_df, totals_all, accepted_all = collect_matches(table, query, gt, affix, rss, extension,
            border, adjust_per_gene_extension, per_gene_prefix, per_gene_suffix)
        summarize_border_stats!(totals_all, accepted_all, extension !== nothing && border > 0)

        isempty(result_df) && return result_df
        raw !== nothing && CSV.write(raw*".gz", result_df, delim='\t', compress=true)

        df = add_counts!(result_df, sequence_lookup)
        sort!(df, [:full_count, :count], rev=[true, true])
        udf = sort(unique(df), [:well, :case, :gene, :db_name, :sequence])
        add_quality_metrics!(udf)

        priority_columns = ["well", "case", "gene", "db_name", "count", "full_count", "ratio", "full_ratio"]
        remaining_columns = setdiff(names(udf), priority_columns)
        udf = udf[:, vcat(priority_columns, remaining_columns)]
        gdf = groupby(udf, [:well, :case, :gene, :db_name, :sequence])
        udf_indexed = transform(gdf, :well => (x -> 1:length(x)) => :flank_index)
        return filter(x -> x.flank_index <= N, udf_indexed)
    end

    # ========================== Reference-gene ratios & lookups ==========================

    function transform_counts(group_df, name; count_col=:count)
        ref_row = filter(row -> startswith(row.db_name, name), group_df)
        ref_count = 1
        well, case = first(map(r->(r.well, r.case), eachrow(unique(group_df, [:well,:case]))))
        # sum() over the reference gene's alleles → a scalar denominator. Using ref_row.count
        # (a vector) errored with DimensionMismatch whenever the refgene matched >1 allele.
        isempty(ref_row) ? (@warn "Reference name $name not found in well $well and case $case") : (@info "Applying name $name to well $well and case $case"; ref_count = sum(ref_row.count))
        group_df[!, "$(count_col)_$(first(split(name,'*')))_ratio"] = group_df[:, count_col] ./ ref_count
        return group_df
    end

    function grouped_ratios(counts_df, refgene; count_col=:count)
        transformed = DataFrame[]
        for group in groupby(counts_df, [:well, :case])
            refgene != "" && (group = transform_counts(group, refgene, count_col=count_col))
            push!(transformed, DataFrame(group))
        end
        return reduce(vcat, transformed)
    end

    function build_sequence_lookup(ref_fasta_path::String)
        sequence_lookup = Dict{String, Bool}()
        @info "Building sequence lookup from reference FASTA: $ref_fasta_path"
        open(FASTA.Reader, ref_fasta_path) do reader
            for record in reader; sequence_lookup[string(FASTA.sequence(record))] = true; end
        end
        @info "Loaded $(length(sequence_lookup)) sequences from reference FASTA"
        return sequence_lookup
    end

    """
        load_ratio_dict(path) -> Dict{String,Float64}

    Load a per-allele/per-gene ratio threshold file (columns `name`, `ratio`) into a typed
    dict. Returns an empty typed dict when `path` is nothing (no throwaway DataFrame).
    """
    function load_ratio_dict(path)
        path === nothing && return Dict{String,Float64}()
        df = CSV.read(path, DataFrame, delim='\t')
        @assert all(n in names(df) for n in ["name", "ratio"]) "ratio file $path must have columns: name, ratio"
        @info "Using ratio file $path with $(nrow(df)) entries"
        return Dict{String,Float64}(string(n) => Float64(r) for (n, r) in zip(df.name, df.ratio))
    end

    # ========================== Output column ordering ==========================
    # Identifiers and metrics (most impactful first) on the left; the long DNA columns (flanks +
    # sequence) on the right in genomic 5'→3' order, so the wide values don't bury the metrics.

    dna_layout(gt::GeneType, ::Integer) = ["prefix", "sequence", "suffix"]   # extension mode
    dna_layout(::VGene, ::Nothing) = ["prefix", "sequence", "heptamer", "spacer", "nonamer"]
    dna_layout(::JGene, ::Nothing) = ["nonamer", "spacer", "heptamer", "sequence", "suffix"]
    dna_layout(::DGene, ::Nothing) = ["pre_nonamer", "pre_spacer", "pre_heptamer", "sequence",
                                      "post_heptamer", "post_spacer", "post_nonamer"]

    const EXACT_LEFT_ORDER = ["well", "case", "gene", "db_name", "isin_db",
        "count", "full_count", "allelic_ratio", "ratio", "full_ratio",
        "n_donors", "n_reads_total", "max_full_ratio",
        "gene_case_freq", "allele_cohort_fold", "allele_cohort_median",
        "gene_cohort_fold", "gene_cohort_median", "gene_count", "case_count",
        "chimera_score", "flank_index", "reject_reason", "reject_stage"]

    """
        order_exact_columns(df, gt, extension) -> df

    Reorder for readability: identifiers + metrics (most impactful first), then any extra columns,
    then the long DNA columns (flanks + sequence) last, in genomic 5'→3' order. Present columns only.
    """
    function order_exact_columns(df::DataFrame, gt::GeneType, extension)
        present = names(df)
        dna = [c for c in dna_layout(gt, extension) if c in present]
        left = [c for c in EXACT_LEFT_ORDER if c in present && !(c in dna)]
        placed = Set(vcat(left, dna))
        middle = [c for c in present if !(c in placed)]
        return select(df, vcat(left, middle, dna))
    end

    # ========================== Findings report ==========================

    "Per-donor (case) depth (reads) and breadth (genes/alleles); flags the weakest donors."
    function report_per_donor(kept::DataFrame, table)
        nrow(kept) == 0 && return nothing
        breadth = combine(groupby(kept, :case),
                          :gene => (x -> length(unique(x))) => :n_genes,
                          :db_name => (x -> length(unique(x))) => :n_alleles)
        # total reads per donor from the demux table; string-keyed so case-id types can differ.
        depth = combine(groupby(table, :case), nrow => :reads)
        reads_by_case = Dict(string(r.case) => r.reads for r in eachrow(depth))
        breadth.reads = [get(reads_by_case, string(c), 0) for c in breadth.case]
        sort!(breadth, :case)
        ndon = nrow(breadth)

        printstyled("  per-donor QC over $ndon donor(s) — genes & reads (low ⇒ donor may have failed):\n";
                    color=:light_black)
        boxplot_if_available(["genes/donor", "alleles/donor"], [breadth.n_genes, breadth.n_alleles])
        boxplot_if_available(["reads/donor"], [breadth.reads])
        worst = first(sort(breadth, :n_genes), min(5, ndon))
        println("    weakest donors by genes detected: ",
                join(["$(r.case): $(r.n_genes)g/$(r.reads)r" for r in eachrow(worst)], ",  "))
        return nothing
    end

    "Per-gene read-count distribution (amplification efficiency), genes sorted by median count."
    function report_per_gene_amplification(kept::DataFrame; max_genes::Int=50)
        nrow(kept) == 0 && return nothing
        genes = String[]; data = Vector{Vector{Int}}()
        for sub in groupby(kept, :gene)
            push!(genes, String(first(sub.gene)))
            push!(data, Vector{Int}(sub.count))
        end
        ord = sortperm([median(d) for d in data]; rev=true)
        genes, data = genes[ord], data[ord]
        if length(genes) > max_genes
            genes, data = genes[1:max_genes], data[1:max_genes]
        end
        printstyled("  per-gene read-count distribution (top = best amplifying; box = spread across alleles/donors):\n";
                    color=:light_black)
        boxplot_if_available(genes, data)
        return nothing
    end

    # Which heptamer column(s) to show, by gene orientation: V's RSS is 3', J's is 5', D has both.
    heptamer_panels(::VGene) = (("heptamer", "heptamer (3' RSS)"),)
    heptamer_panels(::JGene) = (("heptamer", "heptamer (5' RSS)"),)
    heptamer_panels(::DGene) = (("pre_heptamer", "pre-heptamer (5' RSS)"),
                                ("post_heptamer", "post-heptamer (3' RSS)"))

    """
        report_exact_findings(counts_df, kept, db, gt, extension, table)

    Colored diagnostics for an exact search: accepted alleles per gene; per-donor QC (depth &
    breadth, to spot failed donors); per-gene amplification; reference genes never matched or
    fully filtered; the reject-reason breakdown; an accepted-vs-rejected filter-quality comparison;
    and (RSS mode) the heptamer consensus + variation for the side(s) of the searched gene.
    """
    function report_exact_findings(counts_df::DataFrame, kept::DataFrame, db, gt::GeneType, extension, table)
        section("Exact search — findings")
        stage_report("accepted (passed all filters)", nrow(kept), nrow(counts_df))

        query_genes = Set(String(first(split(string(n), '*'))) for (n, _) in db)
        matched = Set(string.(counts_df.gene))
        accepted_g = nrow(kept) > 0 ? Set(string.(kept.gene)) : Set{String}()

        if nrow(kept) > 0
            per_gene = sort(combine(groupby(unique(select(kept, [:gene, :db_name])), :gene),
                                    nrow => :alleles), :alleles, rev=true)
            novel_note = "isin_db" in names(kept) ?
                " ($(count(==("Novel"), string.(kept.isin_db))) novel call row(s))" : ""
            println("  $(sum(per_gene.alleles)) accepted allele(s) across $(nrow(per_gene)) gene(s)$novel_note.")
            printstyled("  accepted alleles per gene:\n"; color=:light_black)
            barplot_if_available(per_gene.gene, per_gene.alleles)
        else
            printstyled("  no alleles passed the filters\n"; color=:light_red)
        end

        report_per_donor(kept, table)
        report_per_gene_amplification(kept)

        never_matched = sort(collect(setdiff(query_genes, matched)))
        fully_filtered = sort(collect(setdiff(matched, accepted_g)))
        if !isempty(never_matched)
            printstyled("  ⚠ $(length(never_matched)) reference gene(s) never matched any read: "; color=:yellow)
            println(join(first(never_matched, 15), ", "), length(never_matched) > 15 ? " …" : "")
        end
        if !isempty(fully_filtered)
            printstyled("  ⚠ $(length(fully_filtered)) gene(s) matched but every candidate was filtered out: "; color=:yellow)
            println(join(first(fully_filtered, 15), ", "), length(fully_filtered) > 15 ? " …" : "")
        end

        report_rejections(counts_df.reject_reason)
        filter_quality_report(counts_df, [:n_donors, :max_full_ratio, :full_count])

        if extension === nothing && nrow(kept) > 0
            for (col, lbl) in heptamer_panels(gt)
                col in names(kept) && rss_consistency(kept[!, Symbol(col)]; label=lbl)
            end
        end
        return nothing
    end

    # ========================== CLI handler (orchestrator) ==========================

    function handle_exact(parsed_args, immunediscover_module, always_gz)
        @info "Exact search"
        ex = parsed_args["search"]["exact"]
        extension = ex["extension"]
        border = get(ex, "border", 0)
        adjust_per_gene_extension = get(ex, "adjust-per-gene-extension", false)
        adjust_percent = get(ex, "adjust-percent", 1.0)
        limit = ex["limit"]
        refgenes = ex["refgene"]
        length(refgenes) > 0 && @info "Using reference genes $refgenes"

        table = immunediscover_module.load_demultiplex(ex["tsv"])
        limit > 0 && (@info "Limiting reads to $limit"; table = table[1:limit, :])
        db = immunediscover_module.load_fasta(ex["fasta"], validate=false)
        mincount = ex["mincount"]
        minratio = ex["minratio"]
        mincount < 5 && @warn "Decreasing mincount below 5 may lead to false positives"
        top = ex["top"]
        affix = ex["affix"]
        gene = ex["gene"]
        locus = ex["locus"]

        local rss
        if extension !== nothing
            @info "Using extension mode with length $extension"; rss = String[]
        else
            rss = split(ex["rss"], ',')
            immunediscover_module.validate_types(rss)
            @info "Extract RSS: $(join(rss,','))"
        end
        top != 1 && @info "Uncollapsed mode; at most $top full records returned."

        # `expect`/`deletion` control-gene threshold files serve two distinct, name-keyed roles:
        # expect_dict also relaxes the within-gene allelic-ratio floor (get_ratio); both feed the
        # allelic_ratio / gene_case_freq floors (get_ratio_threshold). Different thresholds, same list.
        expect_dict = load_ratio_dict(ex["expect"])
        deletion_dict = load_ratio_dict(ex["deletion"])

        raw = ex["raw"]
        sequence_lookup = ex["ref-fasta"] !== nothing ? build_sequence_lookup(ex["ref-fasta"]) : nothing

        counts_df = exact_search(table, db, gene; affix=affix, rss=rss, extension=extension, N=top,
            raw=raw, sequence_lookup=sequence_lookup, border=border,
            adjust_per_gene_extension=adjust_per_gene_extension, adjust_percent=adjust_percent)
        if nrow(counts_df) == 0
            @warn "No exact matches"
            return
        end
        add_chimera_scores!(counts_df, refs_by_gene(db); seq_col=:sequence, gene_col=:gene)
        sort!(counts_df, [:case, :db_name])

        # Count/ratio filters — annotate (don't drop) so the full table records every candidate.
        section("Exact search — count and ratio filters")
        init_rejection_columns!(counts_df)
        annotate_stage!(counts_df,
            exact_filter_criteria(; mincount=mincount, minratio=minratio, expect_dict=expect_dict,
                min_recurrence=get(ex, "min-recurrence", 0), min_seqlen=get(ex, "min-seqlen", 0),
                min_peak_ratio=get(ex, "min-peak-ratio", 0.0)),
            "count and ratio filter")

        if !ex["noplot"]
            plotdf = accepted(counts_df)
            nrow(plotdf) > 0 ? immunediscover_module.plotgenes(plotdf) : @warn "No exact matches to plot"
        end

        @info "Excluding genes not starting with $locus for frequency calculation"
        add_frequency_columns!(counts_df, locus)

        section("Exact search — frequency filters")
        annotate_stage!(counts_df,
            exact_frequency_criteria(immunediscover_module, expect_dict, deletion_dict,
                ex["min-allele-cohort-fold"], ex["min-gene-cohort-fold"]),
            "frequency filter")

        reason_cols = [:reject_reason, :reject_stage]
        kept = accepted(counts_df)
        if length(refgenes) > 0
            for refgene in refgenes
                kept = grouped_ratios(kept, refgene, count_col=:count)
                transform!(groupby(kept, [:well, :case, :gene]), :count => sum => :ref_gene_count)
                kept = grouped_ratios(kept, refgene, count_col=:ref_gene_count)
            end
        end
        output = always_gz(ex["output"])
        full_output = always_gz(replace(replace(output, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv")
        gt = parse_gene_type(gene)
        report_exact_findings(counts_df, kept, db, gt, extension, table)
        # Readable output: metrics left, long sequence/flank columns right (genomic order); round
        # float columns to 4 dp instead of full Float64 precision.
        kept = round_floats!(order_exact_columns(kept, gt, extension))
        counts_df = round_floats!(order_exact_columns(counts_df, gt, extension))
        CSV.write(output, select(kept, Not(reason_cols)), compress=true, delim='\t')
        printstyled("  ✓ "; color=:green, bold=true); println("filtered → $output  ($(nrow(kept)) rows)")
        CSV.write(full_output, counts_df, compress=true, delim='\t')
        printstyled("  ✓ "; color=:green, bold=true); println("full     → $full_output  ($(nrow(counts_df)) candidates + reject_reason)")
        return
    end

    export grouped_ratios, transform_counts, build_sequence_lookup, handle_exact, exact_search
end
