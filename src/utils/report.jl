module Report
    using Statistics
    using Printf
    using ..Align: core_mismatch_row
    using ..Data: histogram_if_available, heatmap_if_available, barplot_if_available

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

    # ===== Reusable findings reporting (rejections, recurrence, filter quality, RSS motif) =====

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
        rss_consistency(seqs; label="heptamer")

    Legible RSS-motif consistency: prints the consensus motif (e.g. `CACAGTG`) and mean
    conservation, then a per-position **variation** bar plot (taller ⇒ that position is less
    conserved), labelled by the consensus base — so it is obvious whether the extracted motif is
    clean and where it varies. A perfectly conserved motif prints a one-line note instead.
    """
    function rss_consistency(seqs; label::AbstractString="heptamer")
        cons, conservation = consensus_motif(seqs)
        isempty(cons) && return nothing
        meanc = mean(conservation)
        printstyled("  ", label, " — consensus ", cons, "  (mean conservation ",
                    round(meanc; digits=3), " = avg fraction of alleles matching the consensus base):\n";
                    color=:light_black)
        variation = round.(1 .- conservation; digits=3)
        if maximum(variation) <= 0.001
            printstyled("    perfectly conserved at every position (variation 0)\n"; color=:light_black)
        else
            printstyled("    per-position variation (0 = fully conserved, taller = more variable; x = position:consensus base):\n";
                        color=:light_black)
            labels = ["$(j):$(cons[j])" for j in 1:length(cons)]
            barplot_if_available(labels, variation)
        end
        return nothing
    end

    const NOVEL_SUFFIX = r"_S\d+$"

    "True when `name` carries a hashed suffix (`IGHV1-2_S1234`)."
    is_novel_name(name) = occursin(NOVEL_SUFFIX, String(name))

    "Hashed discovery name, not an allele already present in the reference DB."
    is_discovered_novel(name, db::AbstractDict{String,String}) =
        is_novel_name(name) && !haskey(db, strip(String(name)))

    "Build a name → sequence lookup from `(name, seq)` DB pairs."
    function db_dict(db_seqs)
        d = Dict{String,String}()
        for (n, s) in db_seqs
            d[strip(String(n))] = String(s)
        end
        return d
    end

    "Full germline sequence for the BLAST best-hit `sseqid` (strips hashed suffix)."
    function matched_germline(sseqid::AbstractString, db::AbstractDict{String,String})
        s = strip(String(sseqid))
        haskey(db, s) && return db[s]
        m = match(r"^(.+)_S\d+$", s)
        m !== nothing && return get(db, String(m.captures[1]), nothing)
        return get(db, s, nothing)
    end

    "SNP weights along `core` vs matched germline (same alignment as `aln_mismatch`)."
    mismatch_row(core, ref, weight) = core_mismatch_row(core, ref, weight)

    """
        snp_support_row(core, ref, nreads, max_reads) -> Vector{Float64}

    Per-position SNP support: 0 at matches, `nreads/max_reads` at mismatches (within-gene
    normalization).
    """
    function snp_support_row(core, ref, nreads::Integer, max_reads::Integer)
        mask = core_mismatch_row(core, ref, 1.0)
        max_reads <= 0 && return zeros(length(mask))
        support = Float64(nreads) / Float64(max_reads)
        return mask .* support
    end

    """
        gene_novel_diff_panels(genes, seqs, names, sseqids; reads, aln_mismatches, db_seqs) ->
            (panels, n_suspicious)

    One matrix per gene with plottable novels. Each row = trimmed core vs its BLAST-matched
    germline allele (`sseqid` in DB). `n_suspicious` counts discovered novels with
    `aln_mismatch == 0` but zero SNP diff on the core (a naming bug).
    """
    function gene_novel_diff_panels(genes, seqs, names, sseqids; reads, aln_mismatches, db_seqs)
        db = db_dict(db_seqs)
        bygene = Dict{String,Dict{String,Tuple{String,Int,String,Int}}}()
        for (g, s, nm, sid, r, mm) in zip(genes, seqs, names, sseqids, reads, aln_mismatches)
            core = String(s)
            isempty(core) && continue
            gene = String(g)
            slot = get!(bygene, gene, Dict{String,Tuple{String,Int,String,Int}}())
            prev = get(slot, core, nothing)
            if prev === nothing || r > prev[2]
                slot[core] = (String(nm), Int(r), String(sid), Int(mm))
            end
        end
        panels = Tuple{String,Matrix{Float64}}[]
        n_suspicious = 0
        for (gene, slot) in bygene
            pending = Tuple{String,String,Int,Int}[]
            for (core, (nm, nreads, sid, mm)) in slot
                is_discovered_novel(nm, db) || continue
                ref = matched_germline(sid, db)
                ref === nothing && continue
                if !any(>(0), mismatch_row(core, ref, 1.0))
                    mm == 0 && (n_suspicious += 1)
                    continue
                end
                push!(pending, (core, ref, nreads, mm))
            end
            isempty(pending) && continue
            max_reads = maximum(p[3] for p in pending)
            row_support(p) = -sum(snp_support_row(p[1], p[2], p[3], max_reads))
            pending = pending[sortperm(pending; by=row_support)]
            rows = [snp_support_row(core, ref, nreads, max_reads) for (core, ref, nreads, _) in pending]
            L = maximum(length(r) for r in rows)
            M = zeros(Float64, length(rows), L)
            for (i, row) in enumerate(rows)
                M[i, 1:length(row)] = row
            end
            push!(panels, (gene, M))
        end
        sort!(panels; by = p -> p[1])
        return panels, n_suspicious
    end

    "Column indices covering every SNP (± `flank` nt), capped for terminal width."
    function snp_window(M::AbstractMatrix{<:Real}; flank::Int=1, max_cols::Int=72)
        cols = Int[]
        for j in axes(M, 2)
            any(>(0), @view(M[:, j])) && append!(cols, max(1, j - flank):min(size(M, 2), j + flank))
        end
        sort!(unique!(cols))
        length(cols) > max_cols && return cols[1:max_cols]
        return cols
    end

    """
        cluster_profile_heatmap(genes, seqs, names, sseqids; reads, parent_ratio, db_seqs, ...)

    Diagnostic heatmap: one panel per gene (sorted alphabetically). Each row is one discovered
    novel trimmed core vs its BLAST-matched germline allele. Color is read-supported SNP signal:
    0 at matches, `n_reads / gene_max_reads` at mismatches (per-position counts are not
    available — support is cluster-level). Likely error satellites are omitted here; see
    `satellite_score`, `likely_satellite`, `nn_dist`, `parent_ratio`, and `chimera_score` in the
    output table. A blank panel means
    every novel in that gene matched germline exactly — investigate naming. Requires `db_seqs`
    and per-row `sseqids`.
    """
    function cluster_profile_heatmap(genes, seqs, names, sseqids; reads, aln_mismatch, db_seqs,
                                     title::AbstractString="novel alleles",
                                     max_rows::Int=12)
        panels, n_suspicious = gene_novel_diff_panels(genes, seqs, names, sseqids;
                                                      reads=reads, aln_mismatches=aln_mismatch,
                                                      db_seqs=db_seqs)
        n_suspicious > 0 && printstyled("      ⚠ ", n_suspicious,
                                        " novel(s) with aln_mismatch=0 but identical to germline — naming bug\n";
                                        color=:yellow)
        isempty(panels) && return nothing
        printstyled("      ", title,
                    " — row = one novel core; color = read-supported SNP ",
                    "(0 = match; brightness ∝ n_reads within gene)\n";
                    color=:light_black)
        for (g, M) in panels
            Mplot = size(M, 1) > max_rows ? M[1:max_rows, :] : M
            L = size(M, 2)
            cols = snp_window(Mplot)
            isempty(cols) && continue
            xlabel = length(cols) < L ? "nt $(first(cols))–$(last(cols)) of $L" : "position"
            nrows = size(M, 1)
            rownote = nrows > max_rows ? ", showing $max_rows/$nrows" : ""
            heatmap_if_available(Mplot[:, cols];
                                 title="$g  ($nrows novel$rownote)",
                                 xlabel=xlabel, ylabel="SNP support", zlim=(0, 1))
        end
        return nothing
    end
end
