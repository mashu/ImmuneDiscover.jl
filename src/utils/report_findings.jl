# Reusable findings: rejections, filter quality, RSS motif consensus.

"""
    reject_counts(reasons; accepted_label="accepted") -> Vector{Tuple{String,Int}}

Tally `reject_reason` strings: the accepted (empty-reason) bucket first, then each reason by
descending count. Pure.
"""
function reject_counts(reasons; accepted_label::AbstractString="accepted")
    d = Dict{String,Int}()
    for r in reasons
        key = isempty(r) ? accepted_label : String(r)
        d[key] = get(d, key, 0) + 1
    end
    rest = sort([(k, v) for (k, v) in d if k != accepted_label]; by = x -> -x[2])
    acc = get(d, accepted_label, 0)
    return acc > 0 ? vcat([(accepted_label, acc)], rest) : rest
end

"Bar plot of how many candidates each filter removed (accepted bucket first)."
function report_rejections(reasons; title::AbstractString="how candidates were filtered")
    counts = reject_counts(reasons)
    isempty(counts) && return nothing
    printstyled("  ", title, " — candidates per outcome:\n"; color=:light_black)
    barplot_if_available([c[1] for c in counts], [c[2] for c in counts])
    return nothing
end

"""
    filter_quality_report(df, metrics; reason_col=:reject_reason)

Compare the median of each metric among ACCEPTED vs REJECTED candidates. When a metric the
filters do not key on (e.g. cross-donor recurrence) is higher for accepted rows, that is
independent evidence the filters keep the stronger candidates. Returns nothing.
"""
function filter_quality_report(df, metrics; reason_col::Symbol=:reject_reason)
    reason_col in propertynames(df) || return nothing
    acc = isempty.(df[!, reason_col])
    (count(acc) == 0 || count(.!acc) == 0) && return nothing
    printstyled("  filter quality — accepted vs rejected medians ",
                "(↑ = filters keep the stronger rows):\n"; color=:light_black)
    for m in metrics
        m in propertynames(df) || continue
        v = df[!, m]
        am = median(v[acc]); rm = median(v[.!acc])
        arrow = am > rm ? "↑" : am < rm ? "↓" : "≈"
        printstyled("    ", rpad(String(m), 16); color=:cyan)
        println("accepted ", round(am; digits=3), "   rejected ", round(rm; digits=3), "   ", arrow)
    end
    return nothing
end

const BASES = ('A', 'C', 'G', 'T')

base_row(::Val{'A'}) = 1
base_row(::Val{'C'}) = 2
base_row(::Val{'G'}) = 3
base_row(::Val{'T'}) = 4
base_row(::Val{C}) where {C} = 0

"""
    dominant_length(seqs) -> Int

Most common sequence length among `seqs` (0 when empty). Used to pick a single, well-defined
alignment width so a base-composition matrix stays rectangular.
"""
function dominant_length(seqs)
    isempty(seqs) && return 0
    counts = Dict{Int,Int}()
    for s in seqs
        counts[length(s)] = get(counts, length(s), 0) + 1
    end
    best_len, best_n = 0, -1
    for (len, n) in counts
        if n > best_n
            best_len, best_n = len, n
        end
    end
    return best_len
end

"""
    composition_matrix(seqs) -> Matrix{Float64}

`4 × L` base-composition matrix (rows A/C/G/T) over equal-length `seqs`; each column is the
per-position base frequency and sums to 1 (ignoring non-ACGT). Conserved columns light up a
single row; variable columns split across rows, localizing where alleles differ. Pure.
"""
function composition_matrix(seqs)
    L = length(first(seqs))
    M = zeros(Float64, 4, L)
    for s in seqs, (j, ch) in enumerate(s)
        i = base_row(Val(ch))
        i == 0 || (M[i, j] += 1.0)
    end
    M ./= length(seqs)
    return M
end

"""
    consensus_motif(seqs) -> (consensus::String, conservation::Vector{Float64})

Per-position consensus base (the most frequent base) and its frequency (conservation, 0–1)
over the dominant-length sequences. Pure.
"""
function consensus_motif(seqs)
    s = [String(x) for x in seqs if !isempty(x)]
    isempty(s) && return ("", Float64[])
    L = dominant_length(s)
    L == 0 && return ("", Float64[])
    keep = [first(x, L) for x in s if length(x) >= L]
    isempty(keep) && return ("", Float64[])
    M = composition_matrix(keep)
    cons = Char[]; conservation = Float64[]
    for j in axes(M, 2)
        col = @view M[:, j]
        i = argmax(col)
        push!(cons, BASES[i]); push!(conservation, col[i])
    end
    return (String(cons), conservation)
end

"""
    rss_consistency(seqs; label="heptamer", color=:green)

Legible RSS-motif consistency: prints the consensus motif (e.g. `CACAGTG`) and mean
conservation, then a per-position **variation** bar plot (taller ⇒ that position is less
conserved), labelled by the consensus base. `color` accents the label and bars so distinct
panels (e.g. a D gene's 5' vs 3' RSS) are visually separable. Perfectly conserved ⇒ one-line note.
"""
function rss_consistency(seqs; label::AbstractString="heptamer", color::Symbol=:green)
    cons, conservation = consensus_motif(seqs)
    isempty(cons) && return nothing
    meanc = mean(conservation)
    printstyled("  ", label; color=color, bold=true)
    printstyled("  — consensus ", cons, "  (mean conservation ", round(meanc; digits=3),
                " = avg fraction of alleles matching the consensus base):\n"; color=:light_black)
    variation = round.(1 .- conservation; digits=3)
    if maximum(variation) <= 0.001
        printstyled("    perfectly conserved at every position (variation 0)\n"; color=:light_black)
    else
        printstyled("    per-position variation (0 = fully conserved, taller = more variable; x = position:consensus base):\n";
                    color=:light_black)
        labels = ["$(j):$(cons[j])" for j in 1:length(cons)]
        barplot_if_available(labels, variation; color=color)
    end
    return nothing
end
