"""
    haplotype

Module for inferring approximate haplotypes from unphased immunoglobulin data.

This module implements haplotype inference based on diploid assumptions, applying
minimum count and allelic ratio thresholds to determine homozygous vs heterozygous
genotypes for each donor and gene. When more than two alleles exceed the allelic
ratio threshold, a duplication genotype is called.

Note: These are approximate haplotypes inferred from unphased data based on 
diploid assumptions. Results should be interpreted with caution.
"""
module haplotype

using DataFrames
using CSV
using Statistics
using FASTX

export infer_haplotypes, handle_haplotype

"""
    load_novel_alleles(fasta_path::String) -> Set{String}

Load novel allele names from a FASTA file.

# Arguments
- `fasta_path`: Path to FASTA file containing novel alleles

# Returns
- Set of novel allele names (sequence headers without '>')
"""
function load_novel_alleles(fasta_path::String)::Set{String}
    novel_alleles = Set{String}()
    
    if isfile(fasta_path)
        reader = FASTA.Reader(open(fasta_path, "r"))
        try
            for record in reader
                allele_name = FASTA.identifier(record)
                push!(novel_alleles, allele_name)
            end
        finally
            close(reader)
        end
    end
    
    return novel_alleles
end


"""
    infer_haplotypes(input_file::String, output_file::String;
                    case_col::String="case", 
                    allele_col::String="db_name",
                    gene_col::String="gene",
                    mincount::Int=5,
                    min_ratio::Float64=0.1,
                    novel_fasta::Union{String, Nothing}=nothing)

Infer approximate haplotypes from unphased immunoglobulin data.

This function groups alleles by gene and donor, applies minimum count and allelic
ratio thresholds, and classifies each donor-gene combination as homozygous or
heterozygous based on diploid assumptions.

# Arguments
- `input_file`: Path to input TSV/TSV.GZ file with allele data
- `output_file`: Path to output TSV file for haplotype results
- `case_col`: Column name for donor/case identifiers
- `allele_col`: Column name for allele identifiers  
- `gene_col`: Column name for gene identifiers
- `mincount`: Minimum count threshold for alleles
- `min_ratio`: Minimum minor/major allele ratio for heterozygous classification
- `novel_fasta`: Optional path to FASTA file with novel alleles

# Output columns
- `case`: Donor identifier
- `gene`: Gene identifier
- `genotype`: "homozygous", "heterozygous", "duplication", or "uncertain"
- `allele_1`: Primary allele (highest count)
- `allele_2`: Secondary allele (second highest count, empty for homozygous)
- `count_1`: Count for primary allele
- `count_2`: Count for secondary allele (0 for homozygous)
- `ratio`: Minor/major allele ratio
- `total_count`: Total count for this gene in this donor
- `other_alleles`: Comma-separated list of additional alleles (when >2 alleles present)
- `novel_1`: Whether allele_1 is novel (if novel_fasta provided)
- `novel_2`: Whether allele_2 is novel (if novel_fasta provided)

# Notes
These are approximate haplotypes inferred from unphased data based on diploid
assumptions. The classification assumes:
- Each donor has at most 2 functional alleles per gene (diploid)
- Allelic ratios below min_ratio indicate homozygous genotype
- Allelic ratios above min_ratio indicate heterozygous genotype
- If three or more alleles (relative to the major) meet min_ratio, we call a
  duplication (multi-copy) event for that gene in the donor
"""
function infer_haplotypes(input_file::String, output_file::String;
                         case_col::String="case", 
                         allele_col::String="db_name",
                         gene_col::String="gene",
                         mincount::Int=5,
                         min_ratio::Float64=0.1,
                         novel_fasta::Union{String, Nothing}=nothing)
    
    # Load input data
    df = CSV.File(input_file, delim='\t') |> DataFrame
    @info "Loaded $(nrow(df)) rows from input file"
    
    # Load novel alleles if provided
    novel_alleles = Set{String}()
    if novel_fasta !== nothing
        novel_alleles = load_novel_alleles(novel_fasta)
        @info "Loaded $(length(novel_alleles)) novel alleles from FASTA"
    end
    
    # Convert column names to symbols
    case_sym = Symbol(case_col)
    allele_sym = Symbol(allele_col)
    gene_sym = Symbol(gene_col)
    
    # Use existing count and ratio data from exact search output
    required_cols = ["count"]
    missing_cols = setdiff(required_cols, names(df))
    if !isempty(missing_cols)
        error("Input file must have columns: $(join(missing_cols, ", ")). This command expects output from exact search or similar commands.")
    end
    
    @info "Using existing count and ratio columns from exact search output"
    
    # Choose count column (prefer full_count if available)
    count_col_sym = :count
    if "full_count" in names(df)
        count_col_sym = :full_count
        @info "Using full_count as allele count source"
    else
        @info "Using count as allele count source"
    end
    
    # Filter by minimum count directly on the input data
    filtered = filter(row -> getproperty(row, count_col_sym) >= mincount, df)
    @info "Filtered to $(nrow(filtered)) allele observations meeting mincount >= $mincount"
    
    # Build extended_sequence to disambiguate variants of the same db_name
    # Use extension mode columns if present; otherwise compose from RSS elements based on gene type
    has_prefix_col = "prefix" in names(filtered)
    has_suffix_col = "suffix" in names(filtered)
    # D gene RSS
    has_pre_nonamer = "pre_nonamer" in names(filtered)
    has_pre_spacer = "pre_spacer" in names(filtered)
    has_pre_heptamer = "pre_heptamer" in names(filtered)
    has_post_heptamer = "post_heptamer" in names(filtered)
    has_post_spacer = "post_spacer" in names(filtered)
    has_post_nonamer = "post_nonamer" in names(filtered)
    # V/J RSS generic names
    has_nonamer = "nonamer" in names(filtered)
    has_spacer = "spacer" in names(filtered)
    has_heptamer = "heptamer" in names(filtered)
    
    function compose_extended_sequence(row)
        gene_id = String(getproperty(row, gene_sym))
        isV = occursin("V", gene_id)
        isD = occursin("D", gene_id)
        isJ = occursin("J", gene_id)
        seq = String(getproperty(row, :sequence))
        pre = ""
        post = ""
        if has_prefix_col
            pre = coalesce(String(getproperty(row, :prefix)), "")
        elseif isJ
            # J gene RSS prefix: nonamer, spacer, heptamer
            if has_nonamer; pre *= coalesce(String(getproperty(row, :nonamer)), ""); end
            if has_spacer; pre *= coalesce(String(getproperty(row, :spacer)), ""); end
            if has_heptamer; pre *= coalesce(String(getproperty(row, :heptamer)), ""); end
        elseif isD
            if has_pre_nonamer; pre *= coalesce(String(getproperty(row, :pre_nonamer)), ""); end
            if has_pre_spacer; pre *= coalesce(String(getproperty(row, :pre_spacer)), ""); end
            if has_pre_heptamer; pre *= coalesce(String(getproperty(row, :pre_heptamer)), ""); end
        end
        if has_suffix_col
            post = coalesce(String(getproperty(row, :suffix)), "")
        elseif isV
            # V gene RSS suffix: heptamer, spacer, nonamer
            if has_heptamer; post *= coalesce(String(getproperty(row, :heptamer)), ""); end
            if has_spacer; post *= coalesce(String(getproperty(row, :spacer)), ""); end
            if has_nonamer; post *= coalesce(String(getproperty(row, :nonamer)), ""); end
        elseif isD
            if has_post_heptamer; post *= coalesce(String(getproperty(row, :post_heptamer)), ""); end
            if has_post_spacer; post *= coalesce(String(getproperty(row, :post_spacer)), ""); end
            if has_post_nonamer; post *= coalesce(String(getproperty(row, :post_nonamer)), ""); end
        end
        return pre * seq * post
    end
    
    filtered[!, :extended_sequence] = [compose_extended_sequence(r) for r in eachrow(filtered)]
    
    # Build stable labels per (db_name, extended_sequence) across all donors by global abundance
    global_agg = combine(groupby(filtered, [allele_sym, :extended_sequence]), count_col_sym => sum => :count_sum)
    label_map = Dict{Tuple{String,String}, String}()
    for g in groupby(global_agg, allele_sym)
        # Sort by global aggregated count desc, tie-break by extended_sequence lexicographically for determinism
        g_sorted = sort(copy(g), [:count_sum, :extended_sequence], rev=[true, false])
        for (rank, row) in enumerate(eachrow(g_sorted))
            base = String(getproperty(row, allele_sym))
            ext = String(row.extended_sequence)
            label_map[(base, ext)] = string(base, "_", rank)
        end
    end
    
    # Group by case and gene to analyze genotypes
    case_gene_groups = groupby(filtered, [case_sym, gene_sym])
    
    results = []
    
    for group_df in case_gene_groups
        case_id = group_df[1, case_sym]
        gene_id = group_df[1, gene_sym]
        
        # Aggregate counts per unique (db_name, extended_sequence) within this case-gene group
        agg = combine(groupby(group_df, [allele_sym, :extended_sequence]), count_col_sym => sum => :count_sum)
        base_names = String.(agg[!, allele_sym])
        exts = String.(agg[!, :extended_sequence])
        counts = agg[!, :count_sum]
        # Stable labels from global mapping
        alleles_display = [get(label_map, (base_names[i], exts[i]), base_names[i]) for i in eachindex(base_names)]
        total_count = sum(counts)
        
        # Sort by aggregated count (descending)
        sorted_indices = sortperm(counts, rev=true)
        sorted_alleles = alleles_display[sorted_indices]
        sorted_counts = counts[sorted_indices]
        
        # Classify genotype based on per-allele ratio to the major allele (from aggregated counts)
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
            selected_alleles = sorted_alleles[1:2]  # keep top 2 in dedicated columns; others go to other_alleles
        else
            genotype_type = "uncertain"
            selected_alleles = sorted_alleles[1:2]
        end
        
        # Prepare output row
        allele_1 = length(selected_alleles) >= 1 ? selected_alleles[1] : ""
        allele_2 = length(selected_alleles) >= 2 ? selected_alleles[2] : ""
        
        # Get counts for selected alleles (already sorted)
        count_1 = length(sorted_counts) >= 1 ? sorted_counts[1] : 0
        count_2 = length(sorted_counts) >= 2 ? sorted_counts[2] : 0
        
        # Calculate ratio
        ratio = count_1 > 0 ? count_2 / count_1 : 0.0
        
        # Capture other alleles
        if genotype_type == "duplication"
            # Only include additional alleles that also meet the min_ratio threshold vs. major
            other_idxs = Int[]
            for i in 3:length(sorted_alleles)
                if ratios_to_major[i] >= min_ratio
                    push!(other_idxs, i)
                end
            end
            other_alleles = isempty(other_idxs) ? "" : join(sorted_alleles[other_idxs], ",")
        else
            # For non-duplication, include any remaining alleles for transparency
            other_alleles = length(sorted_alleles) > 2 ? join(sorted_alleles[3:end], ",") : ""
        end
        
        # Create result tuple with or without novel columns
        freq_1 = total_count > 0 ? count_1 / total_count : 0.0
        freq_2 = total_count > 0 ? count_2 / total_count : 0.0
        if novel_fasta !== nothing
            # Check if alleles are novel using base db_name (independent of suffix)
            base1 = length(sorted_indices) >= 1 ? base_names[sorted_indices[1]] : ""
            base2 = length(sorted_indices) >= 2 ? base_names[sorted_indices[2]] : ""
            novel_1 = base1 in novel_alleles
            novel_2 = base2 in novel_alleles
            
            push!(results, (
                case = case_id,
                gene = gene_id,
                genotype = genotype_type,
                allele_1 = allele_1,
                allele_2 = allele_2,
                count_1 = count_1,
                count_2 = count_2,
                freq_1 = round(freq_1, digits=3),
                freq_2 = round(freq_2, digits=3),
                ratio = round(ratio, digits=3),
                total_count = total_count,
                other_alleles = other_alleles,
                novel_1 = novel_1,
                novel_2 = novel_2
            ))
        else
            push!(results, (
                case = case_id,
                gene = gene_id,
                genotype = genotype_type,
                allele_1 = allele_1,
                allele_2 = allele_2,
                count_1 = count_1,
                count_2 = count_2,
                freq_1 = round(freq_1, digits=3),
                freq_2 = round(freq_2, digits=3),
                ratio = round(ratio, digits=3),
                total_count = total_count,
                other_alleles = other_alleles
            ))
        end
    end
    
    # Create output DataFrame
    results_df = DataFrame(results)
    
    # Save results
    CSV.write(output_file, results_df, delim='\t')
    
    # Summary statistics
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
        
        if novel_fasta !== nothing && heterozygous_count > 0
            novel_in_hetero = sum(skipmissing(results_df[results_df.genotype .== "heterozygous", :novel_1])) + 
                             sum(skipmissing(results_df[results_df.genotype .== "heterozygous", :novel_2]))
            @info "  Novel alleles in heterozygous genotypes: $novel_in_hetero"
        end
    else
        @info "  No genotypes passed the filtering criteria (mincount >= $mincount)"
    end
    
    @info "  Results saved to: $output_file"
    
    return results_df
end

"""
    handle_haplotype(parsed_args)

Handle haplotype command from CLI arguments
"""
function handle_haplotype(parsed_args)
    @info "Haplotype inference"
    input_file = parsed_args["analyze"]["haplotype"]["input"]
    output_file = parsed_args["analyze"]["haplotype"]["output"]
    case_col = parsed_args["analyze"]["haplotype"]["case-col"]
    allele_col = parsed_args["analyze"]["haplotype"]["allele-col"]
    gene_col = parsed_args["analyze"]["haplotype"]["gene-col"]
    mincount = parsed_args["analyze"]["haplotype"]["mincount"]
    min_ratio = parsed_args["analyze"]["haplotype"]["min-ratio"]
    novel_fasta = get(parsed_args["analyze"]["haplotype"], "novel-fasta", nothing)
    
    infer_haplotypes(
        input_file, output_file;
        case_col=case_col,
        allele_col=allele_col,
        gene_col=gene_col,
        mincount=mincount,
        min_ratio=min_ratio,
        novel_fasta=novel_fasta
    )
end

end # module
