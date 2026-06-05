module Report
    using Statistics
    using Printf
    using ..Data: histogram_if_available, heatmap_if_available

    export stage_report, stage_summary, distribution_summary, section, cluster_profile_heatmap

    """
        section(title)

    Print a bold colored section header. Returns nothing.
    """
    function section(title::AbstractString)
        bar = "━"^max(3, 56 - length(title))
        printstyled("\n━━ ", title, " ", bar, "\n"; color=:blue, bold=true)
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
        stage_report(name, kept, before; values=nothing)

    Print a colored per-stage line (respects terminal color). When `values` is given, also
    print a dimmed distribution summary; when `histogram=true`, draw a unicode histogram of
    `values` underneath. Returns nothing.
    """
    function stage_report(name::AbstractString, kept::Integer, before::Integer;
                          values=nothing, histogram::Bool=false)
        removed = before - kept
        printstyled("  ▸ "; color=:magenta, bold=true)
        print(rpad(name, 28), " ")
        printstyled("kept $kept"; color=:cyan, bold=true)
        print("/$before")
        removed > 0 && printstyled("  −$removed"; color=:light_red)
        before > 0 && print("  (", round(100 * kept / before; digits=1), "%)")
        println()
        if values !== nothing && !isempty(values)
            printstyled("      ", distribution_summary(values), "\n"; color=:light_black)
            histogram && histogram_if_available(values; nbins=20, title=name)
        end
        return nothing
    end

    const BASES = ('A', 'C', 'G', 'T')

    "Most common length among the sequences (ties broken arbitrarily)."
    function dominant_length(seqs)
        tally = Dict{Int,Int}()
        for s in seqs
            tally[length(s)] = get(tally, length(s), 0) + 1
        end
        best_len, best_n = 0, -1
        for (len, n) in tally
            n > best_n && (best_len, best_n = len, n)
        end
        return best_len
    end

    """
        cluster_profile_heatmap(seqs; title)

    Draw a 4×L base-composition heatmap (rows A/C/G/T) over the sequences sharing the dominant
    length — a quick view of how consistent / variable a set of candidate sequences is. No-op
    when fewer than two sequences share that length.
    """
    function cluster_profile_heatmap(seqs; title::AbstractString="base composition")
        isempty(seqs) && return nothing
        L = dominant_length(seqs)
        sel = [s for s in seqs if length(s) == L]
        length(sel) < 2 && return nothing
        M = zeros(Float64, 4, L)
        for s in sel, (j, ch) in enumerate(s)
            i = findfirst(==(ch), BASES)
            i === nothing || (M[i, j] += 1.0)
        end
        M ./= length(sel)
        heatmap_if_available(M; title="$title (rows A/C/G/T, n=$(length(sel)), L=$L)")
        return nothing
    end
end
