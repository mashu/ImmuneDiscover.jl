module Report
    using Statistics
    using Printf
    using ..Align: core_mismatch_row
    using ..Data: histogram_if_available, heatmap_if_available, barplot_if_available
    using ..Option: Absent, Present, absent, optional

    export stage_report, stage_summary, distribution_summary, section, cluster_profile_heatmap,
           params_report, reject_counts, report_rejections,
           filter_quality_report, rss_consistency, consensus_motif

    """
        section(title)

    Print a bold colored section header. Returns nothing.
    """
    function section(title::AbstractString)
        bar = "━"^max(3, 56 - length(title))
        printstyled("\n━━ ", title, " ", bar, "\n"; color=:blue, bold=true)
        return nothing
    end

    format_value(::Nothing) = "(unset)"
    format_value(::Absent) = "(unset)"
    format_value(v::AbstractString) = isempty(v) ? "(empty)" : v
    format_value(v) = string(v)

    """
        params_report(block, groups; title)

    Print a CLI argument `block` (a Dict) grouped and ordered by `groups` — a vector of
    `"Group name" => ["key1", "key2", …]` pairs — with a header per group. Any keys not listed
    in `groups` (and not internal `%…%` keys) are printed under "other", so nothing is hidden.
    """
    function params_report(block::AbstractDict, groups; title::AbstractString="parameters")
        section(title)
        shown = Set{String}()
        for (gname, keys) in groups
            present = [(k, block[k]) for k in keys if haskey(block, k)]
            isempty(present) && continue
            printstyled("  ", gname, "\n"; color=:cyan, bold=true)
            for (k, v) in present
                push!(shown, String(k))
                printstyled("    --", rpad(String(k), 18); color=:light_black)
                println(format_value(v))
            end
        end
        leftover = sort([(String(k), v) for (k, v) in block if !(String(k) in shown) && !startswith(String(k), "%")])
        if !isempty(leftover)
            printstyled("  other\n"; color=:cyan, bold=true)
            for (k, v) in leftover
                printstyled("    --", rpad(k, 18); color=:light_black)
                println(format_value(v))
            end
        end
        return nothing
    end

    """
        stage_summary(name, kept, before) -> String

    Plain one-line summary of a pipeline stage (kept/before, dropped, percent). Pure; used by
    `stage_report` and easy to test.
    """
    function stage_summary(name::AbstractString, kept::Integer, before::Integer)
        removed = before - kept
        pct = before > 0 ? round(100 * kept / before; digits=1) : 0.0
        return "$name: kept $kept/$before (−$removed, $pct%)"
    end

    """
        distribution_summary(values) -> String

    One-line min/median/mean/max summary of a numeric collection ("" when empty).
    """
    function distribution_summary(values)
        isempty(values) && return ""
        v = collect(values)
        return @sprintf("min=%.4g  median=%.4g  mean=%.4g  max=%.4g",
                        minimum(v), median(v), mean(v), maximum(v))
    end

    """
        stage_report(name, kept, before; values=(), histogram=false)

    Print a colored per-stage line (respects terminal color). When `values` is non-empty, also
    print a dimmed distribution summary; when `histogram=true`, draw a unicode histogram of
    `values` underneath. Returns nothing.
    """
    function stage_report(name::AbstractString, kept::Integer, before::Integer;
                          values=(), histogram::Bool=false, hist_closed::Symbol=:left)
        removed = before - kept
        printstyled("  ▸ "; color=:magenta, bold=true)
        print(rpad(name, 28), " ")
        printstyled("kept $kept"; color=:cyan, bold=true)
        print("/$before")
        removed > 0 && printstyled("  −$removed"; color=:light_red)
        before > 0 && print("  (", round(100 * kept / before; digits=1), "%)")
        println()
        if !isempty(values)
            printstyled("      ", distribution_summary(values), "\n"; color=:light_black)
            histogram && histogram_if_available(values; nbins=20, title=name, closed=hist_closed)
        end
        return nothing
    end

    include("report_findings.jl")
    include("report_heatmap.jl")
end
