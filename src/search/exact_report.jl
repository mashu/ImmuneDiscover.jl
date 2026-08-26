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
    @info "Building sequence lookup from reference FASTA: $ref_fasta_path"
    records = Data.load_fasta(ref_fasta_path)
    sequence_lookup = Dict(seq => true for (_, seq) in records)
    @info "Loaded $(length(sequence_lookup)) sequences from reference FASTA"
    return sequence_lookup
end

"""
    load_ratio_dict(path) -> Dict{String,Float64}

Load a per-allele/per-gene ratio threshold file (columns `name`, `ratio`) into a typed
dict. Returns an empty typed dict when the path is absent.
"""
load_ratio_dict(path) = load_ratio_dict(optional(path))
load_ratio_dict(::Absent) = Dict{String,Float64}()
function load_ratio_dict(p::Present)
    path = p.value
    df = CSV.read(path, DataFrame, delim='\t')
    @assert all(n in names(df) for n in ["name", "ratio"]) "ratio file $path must have columns: name, ratio"
    @info "Using ratio file $path with $(nrow(df)) entries"
    return Dict{String,Float64}(string(n) => Float64(r) for (n, r) in zip(df.name, df.ratio))
end

# ========================== Output column ordering ==========================
# Identifiers and metrics (most impactful first) on the left; the long DNA columns (flanks +
# sequence) on the right in genomic 5'→3' order, so the wide values don't bury the metrics.

dna_layout(gt::GeneType, ::Integer) = ["prefix", "sequence", "suffix"]   # extension mode
dna_layout(gt::GeneType, e::Present) = dna_layout(gt, e.value)
dna_layout(::VGene, ::Absent) = ["prefix", "sequence", "heptamer", "spacer", "nonamer"]
dna_layout(::JGene, ::Absent) = ["nonamer", "spacer", "heptamer", "sequence", "suffix"]
dna_layout(::DGene, ::Absent) = ["pre_nonamer", "pre_spacer", "pre_heptamer", "sequence",
                                 "post_heptamer", "post_spacer", "post_nonamer"]

const EXACT_LEFT_ORDER = ["well", "case", "gene", "db_name", "isin_db",
    "count", "full_count", "gene_count", "case_count", "n_reads_total", "n_donors",
    "allelic_ratio", "full_allelic_ratio", "peak_allelic_ratio", "gene_fraction",
    "gene_case_freq", "allele_cohort_fold", "gene_cohort_fold",
    "allele_cohort_median", "gene_cohort_median",
    "chimera_score", "flank_index", "reject_reason", "reject_stage"]

# Intermediate statistics omitted from the TSV unless `--diagnostic`. Filters still use them;
# `reject_reason` on the full table is the default way to see why a row dropped.
const EXACT_DIAGNOSTIC_COLUMNS = ["gene_count", "case_count", "n_reads_total", "n_donors",
    "peak_allelic_ratio", "gene_fraction", "gene_case_freq",
    "allele_cohort_fold", "gene_cohort_fold", "allele_cohort_median", "gene_cohort_median",
    "chimera_score", "flank_index", "prefix_len", "suffix_len", "ref_gene_count"]

"""
    order_exact_columns(df, gt, extension) -> df

Reorder for readability: identifiers + metrics (most impactful first), then any extra columns,
then the long DNA columns (flanks + sequence) last, in genomic 5'→3' order. Present columns only.
"""
function order_exact_columns(df::DataFrame, gt::GeneType, extension)
    present = names(df)
    dna = [c for c in dna_layout(gt, optional(extension)) if c in present]
    left = [c for c in EXACT_LEFT_ORDER if c in present && !(c in dna)]
    placed = Set(vcat(left, dna))
    middle = [c for c in present if !(c in placed)]
    return select(df, vcat(left, middle, dna))
end

"""
    select_exact_output_columns(df, gt, extension; diagnostic=false) -> df

`order_exact_columns` then, unless `diagnostic`, drop intermediate statistics
(`EXACT_DIAGNOSTIC_COLUMNS`). `flank_index` stays when `--top` produced more than one
flank variant. Filters still run on the dropped columns.
"""
function select_exact_output_columns(df::DataFrame, gt::GeneType, extension; diagnostic::Bool=false)
    ordered = order_exact_columns(df, gt, extension)
    return drop_exact_diagnostic_columns(ordered, Val(diagnostic))
end

drop_exact_diagnostic_columns(df::DataFrame, ::Val{true}) = df
function drop_exact_diagnostic_columns(df::DataFrame, ::Val{false})
    present = names(df)
    drop = [c for c in EXACT_DIAGNOSTIC_COLUMNS if c in present]
    if "flank_index" in drop && nrow(df) > 0 && maximum(df.flank_index) > 1
        drop = [c for c in drop if c != "flank_index"]
    end
    isempty(drop) && return df
    return select(df, Not(drop))
end

log_exact_column_mode(diagnostic::Bool) = log_exact_column_mode(Val(diagnostic))
log_exact_column_mode(::Val{true}) =
    @info "Writing diagnostic columns (cohort folds, chimera_score, denominators, …)"
log_exact_column_mode(::Val{false}) =
    @info "Slim TSV (identifiers, counts, allelic ratios, sequence/flanks). Pass --diagnostic for intermediate statistics."

add_chimera_if_diagnostic!(df, db, diagnostic::Bool) = add_chimera_if_diagnostic!(df, db, Val(diagnostic))
add_chimera_if_diagnostic!(df, db, ::Val{true}) =
    add_chimera_scores!(df, refs_by_gene(db); seq_col=:sequence, gene_col=:gene)
add_chimera_if_diagnostic!(df, db, ::Val{false}) = df

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

# Which heptamer column(s) to show, by gene orientation, with a side colour (5' = cyan, 3' =
# yellow) so a D gene's two RSS panels are visually separable: V's RSS is 3', J's is 5', D both.
heptamer_panels(::VGene) = (("heptamer", "heptamer (3' RSS)", :yellow),)
heptamer_panels(::JGene) = (("heptamer", "heptamer (5' RSS)", :cyan),)
heptamer_panels(::DGene) = (("pre_heptamer", "pre-heptamer (5' RSS)", :cyan),
                            ("post_heptamer", "post-heptamer (3' RSS)", :yellow))

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
        allele_title = "Accepted alleles per gene"
        printstyled("  ", allele_title, ":\n"; color=:light_black)
        barplot_if_available(per_gene.gene, per_gene.alleles; title=allele_title)
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
    filter_quality_report(counts_df, [:n_donors, PEAK_ALLELIC_RATIO, :full_count])
    report_rss_motifs(optional(extension), kept, gt)
    return nothing
end

report_rss_motifs(::Present, _, _) = nothing
function report_rss_motifs(::Absent, kept, gt)
    nrow(kept) > 0 || return nothing
    for (col, lbl, clr) in heptamer_panels(gt)
        col in names(kept) && rss_consistency(kept[!, Symbol(col)]; label=lbl, color=clr)
    end
    return nothing
end
