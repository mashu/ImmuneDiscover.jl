module Report
    using Statistics
    using Printf

    export stage_report, stage_summary, distribution_summary, section

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

    Print a colored per-stage line (respects terminal color); when `values` is given, also
    print a dimmed distribution summary underneath. Returns nothing.
    """
    function stage_report(name::AbstractString, kept::Integer, before::Integer; values=nothing)
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
        end
        return nothing
    end
end
