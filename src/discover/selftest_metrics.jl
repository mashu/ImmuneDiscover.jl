# Column-schema views: sequence columns contribute their length; others are numeric.
# Functors dispatch on the element, so mixed CSV columns do not need `isa` loops.

struct AsFloat64 end
(::AsFloat64)(::Missing) = NaN
(::AsFloat64)(x::Real) = Float64(x)
(::AsFloat64)(x::AbstractString) = something(tryparse(Float64, String(x)), NaN)
(::AsFloat64)(x) = something(tryparse(Float64, string(x)), NaN)

struct AsSeqLength end
(::AsSeqLength)(::Missing) = 0.0
(::AsSeqLength)(x::AbstractString) = Float64(length(x))
(::AsSeqLength)(x) = Float64(length(string(x)))

metric_values(df::DataFrame, v::LengthView) =
    v.name in names(df) ? AsSeqLength().(df[!, v.name]) : Float64[]
metric_values(df::DataFrame, v::NumericView) =
    v.name in names(df) ? AsFloat64().(df[!, v.name]) : Float64[]

struct KeepGe end
struct KeepLe end
count_keep(::KeepGe, vals, thr) = count(>=(thr), vals)
count_keep(::KeepLe, vals, thr) = count(<=(thr), vals)
dir_label(::KeepGe) = "keep ≥"
dir_label(::KeepLe) = "keep ≤"

function best_threshold(tp::AbstractVector{<:Real}, fp::AbstractVector{<:Real})
    cuts = sort!(unique(vcat(collect(tp), collect(fp))))
    ntp, nfp = length(tp), length(fp)
    best = (KeepGe(), NaN, -Inf, 0, 0)
    for thr in cuts, dir in (KeepGe(), KeepLe())
        keep_tp = count_keep(dir, tp, thr)
        keep_fp = count_keep(dir, fp, thr)
        j = keep_tp / ntp - keep_fp / nfp
        j > best[3] && (best = (dir, Float64(thr), j, keep_tp, keep_fp))
    end
    return best
end

resolved_metrics(discovery, metrics::AbstractVector{<:AbstractString}, blast_block) =
    isempty(metrics) ? blast_discoverable_metrics(discovery, optional(blast_block)) : metrics

function metric_separation(discovery::DataFrame, base::AbstractSet{String},
                           truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                           seq_col::Symbol=:aln_qseq, substring::Bool=true,
                           metrics::AbstractVector{<:AbstractString}=String[],
                           blast_block=absent)
    ms = resolved_metrics(discovery, metrics, blast_block)
    cols = names(discovery)
    String(seq_col) in cols ||
        error("Discovery table has no '$seq_col' column — pass --seq-col, or use the full table (<output>.full.tsv.gz).")
    seqs = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]
    truth_novel = Set(String(s) for (_, s) in truth if is_novel(s, base))

    labels = fill(:none, length(seqs))
    for i in eachindex(seqs)
        c = seqs[i]
        (isempty(c) || !is_novel(c, base)) && continue
        hit = c in truth_novel || (substring && any(t -> seq_match(c, t; substring=true), truth_novel))
        labels[i] = hit ? :tp : :fp
    end
    tp_idx = findall(==(:tp), labels)
    fp_idx = findall(==(:fp), labels)

    Row = NamedTuple{(:metric, :n_tp, :n_fp, :tp_median, :fp_median, :direction,
                      :threshold, :cli_suggestion, :youden, :tp_kept, :fp_removed),
                     Tuple{String,Int,Int,Float64,Float64,String,Float64,String,Float64,Int,Int}}
    rows = Row[]
    for m in ms
        m in cols || continue
        vals = metric_values(discovery, metric_view(m, seq_col))
        tp = filter(isfinite, vals[tp_idx])
        fp = filter(isfinite, vals[fp_idx])
        (isempty(tp) || isempty(fp)) && continue
        dir, thr, j, keep_tp, keep_fp = best_threshold(tp, fp)
        direction = dir_label(dir)
        push!(rows, (metric=m, n_tp=length(tp), n_fp=length(fp),
                     tp_median=median(tp), fp_median=median(fp),
                     direction=direction, threshold=thr,
                     cli_suggestion=blast_cli_suggestion(m, direction, thr),
                     youden=j, tp_kept=keep_tp, fp_removed=length(fp) - keep_fp))
    end
    df = DataFrame(metric=String[], n_tp=Int[], n_fp=Int[], tp_median=Float64[],
                   fp_median=Float64[], direction=String[], threshold=Float64[],
                   cli_suggestion=String[], youden=Float64[], tp_kept=Int[], fp_removed=Int[])
    isempty(rows) || (df = sort!(DataFrame(rows), :youden, rev=true))
    return df
end

function recall_safe_filters(discovery::DataFrame, base::AbstractSet{String},
                             truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                             seq_col::Symbol=:aln_qseq, substring::Bool=true,
                             metrics::AbstractVector{<:AbstractString}=String[],
                             blast_block=absent)
    ms = resolved_metrics(discovery, metrics, blast_block)
    cols = names(discovery)
    String(seq_col) in cols ||
        error("Discovery table has no '$seq_col' column — pass --seq-col, or use the full table (<output>.full.tsv.gz).")
    getcol(name) = name in cols ? [ismissing(x) ? "" : String(x) for x in discovery[!, name]] :
                                  fill("", nrow(discovery))
    reason = getcol("reject_reason")
    seqs   = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]
    acc = isempty.(reason)

    truth_novel = [(String(n), String(s)) for (n, s) in truth if is_novel(s, base)]

    core_rows = Dict{String,Vector{Int}}()
    for i in eachindex(seqs)
        (acc[i] && !isempty(seqs[i]) && is_novel(seqs[i], base)) || continue
        push!(get!(core_rows, seqs[i], Int[]), i)
    end
    isempty(core_rows) && return empty_safe_df()
    cores = collect(keys(core_rows))
    is_tp = Dict(c => any(t -> seq_match(c, t[2]; substring=substring), truth_novel) for c in cores)
    allele_cores = Dict{String,Vector{String}}()
    for (name, s) in truth_novel
        matched = [c for c in cores if seq_match(c, s; substring=substring)]
        isempty(matched) || (allele_cores[name] = matched)
    end
    fp_cores = [c for c in cores if !is_tp[c]]

    Row = NamedTuple{(:metric, :acc_fp, :direction, :threshold, :cli_suggestion, :fp_removed),
                     Tuple{String,Int,String,Float64,String,Int}}
    rows = Row[]
    for m in ms
        m in cols || continue
        vals = metric_values(discovery, metric_view(m, seq_col))
        cmax = Dict{String,Float64}(); cmin = Dict{String,Float64}()
        for (c, idx) in core_rows
            fv = filter(isfinite, vals[idx])
            isempty(fv) && continue
            cmax[c] = maximum(fv); cmin[c] = minimum(fv)
        end
        (isempty(allele_cores) || isempty(fp_cores)) && continue
        allele_ge = [maximum(cmax[c] for c in cs if haskey(cmax, c); init=-Inf) for cs in values(allele_cores)]
        allele_le = [minimum(cmin[c] for c in cs if haskey(cmin, c); init=Inf)  for cs in values(allele_cores)]
        safe_ge = minimum(allele_ge); safe_le = maximum(allele_le)
        fp_ge = count(c -> haskey(cmax, c) && cmax[c] < safe_ge, fp_cores)
        fp_le = count(c -> haskey(cmin, c) && cmin[c] > safe_le, fp_cores)
        if fp_ge >= fp_le
            direction = dir_label(KeepGe())
            push!(rows, (metric=m, acc_fp=length(fp_cores), direction=direction, threshold=safe_ge,
                         cli_suggestion=blast_cli_suggestion(m, direction, safe_ge),
                         fp_removed=fp_ge))
        else
            direction = dir_label(KeepLe())
            push!(rows, (metric=m, acc_fp=length(fp_cores), direction=direction, threshold=safe_le,
                         cli_suggestion=blast_cli_suggestion(m, direction, safe_le),
                         fp_removed=fp_le))
        end
    end
    isempty(rows) && return empty_safe_df()
    return sort!(DataFrame(rows), :fp_removed, rev=true)
end

empty_safe_df() = DataFrame(metric=String[], acc_fp=Int[], direction=String[],
                             threshold=Float64[], cli_suggestion=String[], fp_removed=Int[])

function marginal_filter_shadowing(discovery::DataFrame, base::AbstractSet{String},
                                  truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                                  seq_col::Symbol=:aln_qseq, substring::Bool=true,
                                  blast_block=absent)
    return shadow_filters(optional(blast_block), discovery, base, truth, seq_col, substring)
end

shadow_filters(::Absent, _, _, _, _, _) = Dict{String,Int}()
shadow_filters(b::Present, discovery, base, truth, seq_col, substring) =
    shadow_filters(b.value, discovery, base, truth, seq_col, substring)

function shadow_filters(blast_block::AbstractDict, discovery, base, truth, seq_col, substring)
    cols = names(discovery)
    String(seq_col) in cols || return Dict{String,Int}()
    getcol(name) = name in cols ? [ismissing(x) ? "" : String(x) for x in discovery[!, name]] :
                                  fill("", nrow(discovery))
    reason = getcol("reject_reason")
    seqs = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]
    criteria = build_blast_output_criteria(blast_block)
    truth_novel = [(String(n), String(s)) for (n, s) in truth if is_novel(s, base)]
    shadow = Dict{String,Int}()
    for (name, s) in truth_novel
        for i in eachindex(seqs)
            isempty(reason[i]) && continue
            seq_match(seqs[i], s; substring=substring) || continue
            recorded = reason[i]
            row = discovery[i, :]
            for c in criteria
                passes(row, c) && continue
                c.label == recorded && continue
                shadow[c.label] = get(shadow, c.label, 0) + 1
            end
        end
    end
    return shadow
end
