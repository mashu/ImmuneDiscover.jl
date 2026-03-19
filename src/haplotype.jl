"""
    haplotype

Module for inferring approximate haplotypes from unphased immunoglobulin data.
"""
module Haplotype

using DataFrames
using CSV
using Statistics
using FASTX
using ..Exact: GeneType, VGene, DGene, JGene, gene_type_from_name

export infer_haplotypes, handle_haplotype

const HaplotypeRow = NamedTuple{
    (:case,:gene,:genotype,:allele_1,:allele_2,:count_1,:count_2,:freq_1,:freq_2,:ratio,:total_count,:other_alleles,:novel_1,:novel_2),
    Tuple{String,String,String,String,String,Int,Int,Float64,Float64,Float64,Int,String,Bool,Bool}}

function load_novel_alleles(fasta_path::String)::Set{String}
    novel_alleles = Set{String}()
    if isfile(fasta_path)
        open(FASTA.Reader, fasta_path) do reader
            for record in reader
                push!(novel_alleles, FASTA.identifier(record))
            end
        end
    end
    return novel_alleles
end

# Handle missing values safely — String(missing) throws before coalesce sees it.
function safe_string(row, col::Symbol)
    val = getproperty(row, col)
    ismissing(val) ? "" : String(val)
end

# --- Flank collection via dispatch (replaces occursin("V/D/J", gene_id) conditionals) ---

"""
Tracks which flank columns are available in the DataFrame.
Avoids repeated `in names(df)` checks per row.
"""
struct FlankColumns
    prefix::Bool
    suffix::Bool
    heptamer::Bool
    spacer::Bool
    nonamer::Bool
    pre_heptamer::Bool
    pre_spacer::Bool
    pre_nonamer::Bool
    post_heptamer::Bool
    post_spacer::Bool
    post_nonamer::Bool
end

function FlankColumns(df::DataFrame)
    FlankColumns(
        "prefix" in names(df), "suffix" in names(df),
        "heptamer" in names(df), "spacer" in names(df), "nonamer" in names(df),
        "pre_heptamer" in names(df), "pre_spacer" in names(df), "pre_nonamer" in names(df),
        "post_heptamer" in names(df), "post_spacer" in names(df), "post_nonamer" in names(df))
end

# Default: prefix column if available, otherwise empty
function collect_pre_flank(row, fc::FlankColumns, ::GeneType)
    fc.prefix ? safe_string(row, :prefix) : ""
end

# Default: suffix column if available, otherwise empty
function collect_post_flank(row, fc::FlankColumns, ::GeneType)
    fc.suffix ? safe_string(row, :suffix) : ""
end

# V genes: post-flank is heptamer+spacer+nonamer (RSS after V segment)
function collect_post_flank(row, fc::FlankColumns, ::VGene)
    fc.suffix && return safe_string(row, :suffix)
    buf = IOBuffer()
    fc.heptamer && print(buf, safe_string(row, :heptamer))
    fc.spacer && print(buf, safe_string(row, :spacer))
    fc.nonamer && print(buf, safe_string(row, :nonamer))
    String(take!(buf))
end

# J genes: pre-flank is nonamer+spacer+heptamer (RSS before J segment)
function collect_pre_flank(row, fc::FlankColumns, ::JGene)
    fc.prefix && return safe_string(row, :prefix)
    buf = IOBuffer()
    fc.nonamer && print(buf, safe_string(row, :nonamer))
    fc.spacer && print(buf, safe_string(row, :spacer))
    fc.heptamer && print(buf, safe_string(row, :heptamer))
    String(take!(buf))
end

# D genes: pre-flank is pre_nonamer+pre_spacer+pre_heptamer
function collect_pre_flank(row, fc::FlankColumns, ::DGene)
    fc.prefix && return safe_string(row, :prefix)
    buf = IOBuffer()
    fc.pre_nonamer && print(buf, safe_string(row, :pre_nonamer))
    fc.pre_spacer && print(buf, safe_string(row, :pre_spacer))
    fc.pre_heptamer && print(buf, safe_string(row, :pre_heptamer))
    String(take!(buf))
end

# D genes: post-flank is post_heptamer+post_spacer+post_nonamer
function collect_post_flank(row, fc::FlankColumns, ::DGene)
    fc.suffix && return safe_string(row, :suffix)
    buf = IOBuffer()
    fc.post_heptamer && print(buf, safe_string(row, :post_heptamer))
    fc.post_spacer && print(buf, safe_string(row, :post_spacer))
    fc.post_nonamer && print(buf, safe_string(row, :post_nonamer))
    String(take!(buf))
end

function compose_extended_sequence(row, gene_sym::Symbol, fc::FlankColumns)
    gene_id = safe_string(row, gene_sym)
    gt = gene_type_from_name(gene_id)
    seq = safe_string(row, :sequence)
    pre = gt !== nothing ? collect_pre_flank(row, fc, gt) : ""
    post = gt !== nothing ? collect_post_flank(row, fc, gt) : ""
    return pre * seq * post
end

function infer_haplotypes(input_file::String, output_file::String;
                         case_col::String="case",
                         allele_col::String="db_name",
                         gene_col::String="gene",
                         mincount::Int=5,
                         min_ratio::Float64=0.1,
                         novel_fasta::Union{String, Nothing}=nothing)

    df = CSV.File(input_file, delim='\t') |> DataFrame
    @info "Loaded $(nrow(df)) rows from input file"

    novel_alleles = Set{String}()
    if novel_fasta !== nothing
        novel_alleles = load_novel_alleles(novel_fasta)
        @info "Loaded $(length(novel_alleles)) novel alleles from FASTA"
    end

    case_sym = Symbol(case_col)
    allele_sym = Symbol(allele_col)
    gene_sym = Symbol(gene_col)

    required_cols = ["count"]
    missing_cols = setdiff(required_cols, names(df))
    if !isempty(missing_cols)
        error("Input file must have columns: $(join(missing_cols, ", ")). This command expects output from exact search or similar commands.")
    end

    @info "Using existing count and ratio columns from exact search output"

    count_col_sym = :count
    if "full_count" in names(df)
        count_col_sym = :full_count
        @info "Using full_count as allele count source"
    else
        @info "Using count as allele count source"
    end

    filtered = filter(row -> getproperty(row, count_col_sym) >= mincount, df)
    @info "Filtered to $(nrow(filtered)) allele observations meeting mincount >= $mincount"

    # Build extended_sequence to disambiguate variants of the same db_name
    fc = FlankColumns(filtered)

    filtered[!, :extended_sequence] = [compose_extended_sequence(r, gene_sym, fc) for r in eachrow(filtered)]

    # Build stable labels per (db_name, extended_sequence) across all donors
    global_agg = combine(groupby(filtered, [allele_sym, :extended_sequence]), count_col_sym => sum => :count_sum)
    label_map = Dict{Tuple{String,String}, String}()
    for g in groupby(global_agg, allele_sym)
        g_sorted = sort(copy(g), [:count_sum, :extended_sequence], rev=[true, false])
        for (rank, row) in enumerate(eachrow(g_sorted))
            base = String(getproperty(row, allele_sym))
            ext = String(row.extended_sequence)
            label_map[(base, ext)] = string(base, "_", rank)
        end
    end

    case_gene_groups = groupby(filtered, [case_sym, gene_sym])

    # Typed result accumulator — consistent shape with novel fields always present
    results = HaplotypeRow[]

    for group_df in case_gene_groups
        case_id = group_df[1, case_sym]
        gene_id = group_df[1, gene_sym]

        agg = combine(groupby(group_df, [allele_sym, :extended_sequence]), count_col_sym => sum => :count_sum)
        base_names = String.(agg[!, allele_sym])
        exts = String.(agg[!, :extended_sequence])
        counts = agg[!, :count_sum]
        alleles_display = [get(label_map, (base_names[i], exts[i]), base_names[i]) for i in eachindex(base_names)]
        total_count = sum(counts)

        sorted_indices = sortperm(counts, rev=true)
        sorted_alleles = alleles_display[sorted_indices]
        sorted_counts = counts[sorted_indices]

        first_count = sorted_counts[1]
        ratios_to_major = first_count > 0 ? sorted_counts ./ first_count : zeros(length(sorted_counts))
        effective = count(x -> x >= min_ratio, ratios_to_major)

        if effective <= 1
            genotype_type = "homozygous"
            selected_alleles = [sorted_alleles[1]]
        elseif effective == 2
            genotype_type = "heterozygous"
            selected_alleles = sorted_alleles[1:2]
        elseif effective >= 3
            genotype_type = "duplication"
            selected_alleles = sorted_alleles[1:2]
        else
            genotype_type = "uncertain"
            selected_alleles = sorted_alleles[1:min(2, length(sorted_alleles))]
        end

        allele_1 = length(selected_alleles) >= 1 ? selected_alleles[1] : ""
        allele_2 = length(selected_alleles) >= 2 ? selected_alleles[2] : ""
        count_1 = length(sorted_counts) >= 1 ? sorted_counts[1] : 0
        count_2 = length(sorted_counts) >= 2 ? sorted_counts[2] : 0
        ratio = count_1 > 0 ? count_2 / count_1 : 0.0

        if genotype_type == "duplication"
            other_idxs = [i for i in 3:length(sorted_alleles) if ratios_to_major[i] >= min_ratio]
            other_alleles = isempty(other_idxs) ? "" : join(sorted_alleles[other_idxs], ",")
        else
            other_alleles = length(sorted_alleles) > 2 ? join(sorted_alleles[3:end], ",") : ""
        end

        freq_1 = total_count > 0 ? count_1 / total_count : 0.0
        freq_2 = total_count > 0 ? count_2 / total_count : 0.0

        # Always include novel_1/novel_2 for consistent NamedTuple shape (type stability)
        base1 = length(sorted_indices) >= 1 ? base_names[sorted_indices[1]] : ""
        base2 = length(sorted_indices) >= 2 ? base_names[sorted_indices[2]] : ""
        is_novel_1 = novel_fasta !== nothing ? (base1 in novel_alleles) : false
        is_novel_2 = novel_fasta !== nothing ? (base2 in novel_alleles) : false

        push!(results, (case=String(case_id), gene=String(gene_id), genotype=genotype_type,
                        allele_1=allele_1, allele_2=allele_2,
                        count_1=count_1, count_2=count_2,
                        freq_1=round(freq_1, digits=3), freq_2=round(freq_2, digits=3),
                        ratio=round(ratio, digits=3), total_count=total_count,
                        other_alleles=other_alleles,
                        novel_1=is_novel_1, novel_2=is_novel_2))
    end

    results_df = DataFrame(results)
    CSV.write(output_file, results_df, delim='\t')

    total_genotypes = nrow(results_df)
    @info "Haplotype inference completed:"
    @info "  Total case-gene combinations: $total_genotypes"

    if total_genotypes > 0
        homozygous_count = count(x -> x == "homozygous", results_df.genotype)
        heterozygous_count = count(x -> x == "heterozygous", results_df.genotype)
        duplication_count = count(x -> x == "duplication", results_df.genotype)
        uncertain_count = count(x -> x == "uncertain", results_df.genotype)
        @info "  Homozygous: $homozygous_count ($(round(100*homozygous_count/total_genotypes, digits=1))%)"
        @info "  Heterozygous: $heterozygous_count ($(round(100*heterozygous_count/total_genotypes, digits=1))%)"
        @info "  Duplication: $duplication_count ($(round(100*duplication_count/total_genotypes, digits=1))%)"
        @info "  Uncertain: $uncertain_count ($(round(100*uncertain_count/total_genotypes, digits=1))%)"
    end
    @info "  Results saved to: $output_file"
    return results_df
end

function handle_haplotype(parsed_args)
    @info "Haplotype inference"
    block = parsed_args["analyze"]["haplotype"]
    infer_haplotypes(block["input"], block["output"];
        case_col=block["case-col"], allele_col=block["allele-col"],
        gene_col=block["gene-col"], mincount=block["mincount"],
        min_ratio=block["min-ratio"],
        novel_fasta=get(block, "novel-fasta", nothing))
end

end # module
