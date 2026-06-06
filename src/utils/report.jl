module Report
    using Statistics
    using Printf
    using ..Data: histogram_if_available, heatmap_if_available

    export stage_report, stage_summary, distribution_summary, section, cluster_profile_heatmap,
           params_report

    """
        section(title)

    Print a bold colored section header. Returns nothing.
    """
    function section(title::AbstractString)
        bar = "━"^max(3, 56 - length(title))
        printstyled("\n━━ ", title, " ", bar, "\n"; color=:blue, bold=true)
        return nothing
    end

    format_value(v) = v === nothing ? "(unset)" : (v == "" ? "(empty)" : string(v))

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
        stage_report(name, kept, before; values=nothing)

    Print a colored per-stage line (respects terminal color). When `values` is given, also
    print a dimmed distribution summary; when `histogram=true`, draw a unicode histogram of
    `values` underneath. Returns nothing.
    """
    function stage_report(name::AbstractString, kept::Integer, before::Integer;
                          values=nothing, histogram::Bool=false, hist_closed::Symbol=:left)
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
            histogram && histogram_if_available(values; nbins=20, title=name, closed=hist_closed)
        end
        return nothing
    end

    const BASES = ('A', 'C', 'G', 'T')

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
            i = findfirst(==(ch), BASES)
            i === nothing || (M[i, j] += 1.0)
        end
        M ./= length(seqs)
        return M
    end

    const NOVEL_SUFFIX = r"_S\d+$"

    "True when `name` is a hashed novel call (`IGHV1-2_S1234`), not a known allele label."
    is_novel_name(name) = occursin(NOVEL_SUFFIX, String(name))

    "Hamming distance for equal-length cores; Levenshtein otherwise."
    function core_distance(a::AbstractString, b::AbstractString)
        la, lb = length(a), length(b)
        if la == lb
            return sum(ca != cb for (ca, cb) in zip(a, b); init=0)
        end
        prev = collect(0:lb)
        cur = similar(prev)
        for i in 1:la
            cur[1] = i
            for j in 1:lb
                cost = a[i] == b[j] ? 0 : 1
                cur[j + 1] = min(prev[j + 1] + 1, cur[j] + 1, prev[j] + cost)
            end
            prev, cur = cur, prev
        end
        return prev[lb + 1]
    end

    """
        neighbor_parents(cores, reads; max_parents) -> Vector{Union{Nothing,String}}

    For each core, the nearest strictly more-abundant core in the same gene (its likely parent).
    Mirrors the BLAST `neighbor_stats` logic; `nothing` when no parent exists.
    """
    function neighbor_parents(cores::AbstractVector{<:AbstractString},
                              reads::AbstractVector{<:Integer};
                              max_parents::Int=128)
        n = length(cores)
        parents = Vector{Union{Nothing,String}}(nothing, n)
        order = sortperm(reads, rev=true)
        for ii in 1:n
            i = order[ii]
            best_d = typemax(Int)
            best_j = 0
            considered = 0
            for jj in 1:(ii - 1)
                j = order[jj]
                reads[j] > reads[i] || continue
                considered += 1
                considered > max_parents && break
                d = core_distance(cores[i], cores[j])
                if d < best_d || (d == best_d && reads[j] > reads[best_j])
                    best_d, best_j = d, j
                end
            end
            best_j > 0 && (parents[i] = String(cores[best_j]))
        end
        return parents
    end

    "Row brightness from read support; dimmer when the core is a likely error satellite."
    function row_confidence(reads::Real, parent_ratio::Real)
        read_w = clamp(log10(reads + 1) / 4.0, 0.1, 1.0)
        sat_w = parent_ratio >= 20 ? 0.2 : parent_ratio >= 5 ? 0.5 : 1.0
        return read_w * sat_w
    end

    "Per-position mismatch weight vs equal-length `ref` (0 = match, `weight` = SNP)."
    function mismatch_row(seq::AbstractString, ref::AbstractString, weight::Real)
        length(seq) == length(ref) || return zeros(Float64, 0)
        row = zeros(Float64, length(seq))
        w = Float64(weight)
        for j in eachindex(row)
            seq[j] != ref[j] && (row[j] = w)
        end
        return row
    end

    "Reference core for diffing: same-length parent, else same-length known, else same-length abundant core."
    function reference_core(core, cores, names_g, reads_by_core, parent_of)
        same_len(c) = length(c) == length(core)
        parent = get(parent_of, core, nothing)
        parent !== nothing && same_len(parent) && return parent
        known = [c for c in cores if c != core && same_len(c) && !is_novel_name(names_g[c])]
        !isempty(known) && return known[argmax(reads_by_core[c] for c in known)]
        others = [c for c in cores if c != core && same_len(c)]
        isempty(others) && return nothing
        return others[argmax(reads_by_core[c] for c in others)]
    end

    """
        gene_novel_diff_panels(genes, seqs, names; reads, parent_ratios) ->
            Vector{Tuple{String,Matrix{Float64}}}

    One `n×L` matrix per gene with ≥1 plottable novel: row = novel vs same-length parent,
    cell = SNP weight (0 = match, brighter = mismatch; dimmer rows = low reads / satellites).
    """
    function gene_novel_diff_panels(genes, seqs, names; reads, parent_ratios)
        bygene = Dict{String,Dict{String,Tuple{String,Int,Float64}}}()
        for (g, s, nm, r, pr) in zip(genes, seqs, names, reads, parent_ratios)
            core = String(s)
            isempty(core) && continue
            gene = String(g)
            slot = get!(bygene, gene, Dict{String,Tuple{String,Int,Float64}}())
            prev = get(slot, core, nothing)
            if prev === nothing || r > prev[2]
                slot[core] = (String(nm), Int(r), Float64(pr))
            end
        end
        panels = Tuple{String,Matrix{Float64}}[]
        for (gene, slot) in bygene
            cores = collect(keys(slot))
            reads_g = [slot[c][2] for c in cores]
            names_g = Dict(c => slot[c][1] for c in cores)
            pr_g = Dict(c => slot[c][3] for c in cores)
            parents = neighbor_parents(cores, reads_g)
            parent_of = Dict(cores[i] => parents[i] for i in eachindex(cores))
            novel_cores = [c for c in cores if is_novel_name(names_g[c])]
            isempty(novel_cores) && continue
            reads_by_core = Dict(c => slot[c][2] for c in cores)
            rows = Tuple{Vector{Float64},Bool,Int}[]
            for core in novel_cores
                ref = reference_core(core, cores, names_g, reads_by_core, parent_of)
                ref === nothing && continue
                length(core) == length(ref) || continue
                row = mismatch_row(core, ref, row_confidence(reads_by_core[core], pr_g[core]))
                isempty(row) && continue
                push!(rows, (row, pr_g[core] >= 20, reads_by_core[core]))
            end
            isempty(rows) && continue
            sort!(rows; by = r -> (r[2], -r[3]))   # satellites last, then by reads
            L = maximum(length(r[1]) for r in rows)
            M = zeros(Float64, length(rows), L)
            for (i, (row, _, _)) in enumerate(rows)
                M[i, 1:length(row)] = row
            end
            push!(panels, (gene, M))
        end
        sort!(panels; by = p -> -sum(p[2]))
        return panels
    end

    """
        cluster_profile_heatmap(genes, seqs; names, reads, parent_ratio, title, max_genes)

    Quick diagnostic: one heatmap per gene with novel candidates. Each row is a novel vs its
    parent core; dark = conserved, bright = SNP; dim rows = satellites / low support. No text
    listing of positions — the plot is the summary.
    """
    function cluster_profile_heatmap(genes, seqs; names, reads, parent_ratio,
                                     title::AbstractString="novel alleles",
                                     max_genes::Int=15, max_rows::Int=12)
        panels = gene_novel_diff_panels(genes, seqs, names; reads=reads, parent_ratios=parent_ratio)
        isempty(panels) && return nothing
        shown = first(panels, min(length(panels), max_genes))
        printstyled("      ", title,
                    " — row = novel vs parent; dark = match, bright = SNP (dim = satellite)\n";
                    color=:light_black)
        for (g, M) in shown
            Mplot = size(M, 1) > max_rows ? M[1:max_rows, :] : M
            heatmap_if_available(Mplot;
                                 title="$g  ($(size(M, 1)) novel, L=$(size(M, 2)) nt)",
                                 xlabel="position", ylabel="")
        end
        extra = length(panels) - length(shown)
        extra > 0 && printstyled("      … ", extra, " more gene(s) not shown\n"; color=:light_black)
        return nothing
    end
end
