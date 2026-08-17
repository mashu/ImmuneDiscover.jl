# Novel-vs-germline SNP heatmaps (diagnostic panels for discover blast).

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
    haskey(db, s) && return Present(db[s])
    m = match(r"^(.+)_S\d+$", s)
    m === nothing && return optional(get(db, s, nothing))
    return optional(get(db, String(m.captures[1]), nothing))
end

"SNP weights along `core` vs matched germline (same alignment as `core_aln_mismatch`)."
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

take_better_core(::Absent, nm, r, sid, mm) = (String(nm), Int(r), String(sid), Int(mm))
function take_better_core(prev::Present, nm, r, sid, mm)
    r > prev.value[2] ? (String(nm), Int(r), String(sid), Int(mm)) : prev.value
end

"""
    gene_novel_diff_panels(genes, seqs, names, sseqids; reads, core_aln_mismatches, db_seqs) ->
        (panels, n_suspicious)

One matrix per gene with plottable novels. Each row = trimmed core vs its BLAST-matched
germline allele (`sseqid` in DB). `n_suspicious` counts discovered novels with
`core_aln_mismatch == 0` but zero SNP diff on the core (a naming bug).
"""
function gene_novel_diff_panels(genes, seqs, names, sseqids; reads, core_aln_mismatches, db_seqs)
    db = db_dict(db_seqs)
    bygene = Dict{String,Dict{String,Tuple{String,Int,String,Int}}}()
    for (g, s, nm, sid, r, mm) in zip(genes, seqs, names, sseqids, reads, core_aln_mismatches)
        core = String(s)
        isempty(core) && continue
        gene = String(g)
        slot = get!(bygene, gene, Dict{String,Tuple{String,Int,String,Int}}())
        prev = optional(get(slot, core, nothing))
        slot[core] = take_better_core(prev, nm, r, sid, mm)
    end
    panels = Tuple{String,Matrix{Float64}}[]
    n_suspicious = 0
    for (gene, slot) in bygene
        pending = Tuple{String,String,Int,Int}[]
        n_suspicious += collect_gene_novels!(pending, slot, db)
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

collect_gene_novels!(pending, slot, db) = begin
    n = 0
    for (core, rec) in slot
        n += maybe_push_novel!(pending, core, rec, db)
    end
    return n
end

function maybe_push_novel!(pending, core, rec, db)
    nm, nreads, sid, mm = rec
    is_discovered_novel(nm, db) || return 0
    return push_if_ref(pending, core, nreads, mm, matched_germline(sid, db))
end

push_if_ref(_, _, _, _, ::Absent) = 0
function push_if_ref(pending, core, nreads, mm, ref::Present)
    if !any(>(0), mismatch_row(core, ref.value, 1.0))
        return Int(mm == 0)
    end
    push!(pending, (core, ref.value, nreads, mm))
    return 0
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
function cluster_profile_heatmap(genes, seqs, names, sseqids; reads, core_aln_mismatch, db_seqs,
                                 title::AbstractString="novel alleles",
                                 max_rows::Int=12)
    panels, n_suspicious = gene_novel_diff_panels(genes, seqs, names, sseqids;
                                                  reads=reads,
                                                  core_aln_mismatches=core_aln_mismatch,
                                                  db_seqs=db_seqs)
    n_suspicious > 0 && printstyled("      ⚠ ", n_suspicious,
                                    " novel(s) with core_aln_mismatch=0 but identical to germline — naming bug\n";
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
