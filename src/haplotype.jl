"""
    haplotype

Module for inferring approximate haplotypes from unphased immunoglobulin data.

This module implements haplotype inference based on diploid assumptions, applying
minimum count and allelic ratio thresholds to determine homozygous vs heterozygous
genotypes for each donor and gene.

Note: These are approximate haplotypes inferred from unphased data based on 
diploid assumptions. Results should be interpreted with caution.
"""
module haplotype

using DataFrames
using CSV
using Statistics
using FASTX

export infer_haplotypes

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
- `genotype`: "homozygous", "heterozygous", or "uncertain"
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
    
    # Filter by minimum count directly on the input data
    filtered = filter(row -> row.count >= mincount, df)
    @info "Filtered to $(nrow(filtered)) allele observations meeting mincount >= $mincount"
    
    # Group by case and gene to analyze genotypes
    case_gene_groups = groupby(filtered, [case_sym, gene_sym])
    
    results = []
    
    for group_df in case_gene_groups
        case_id = group_df[1, case_sym]
        gene_id = group_df[1, gene_sym]
        
        # Extract alleles and counts for this case-gene combination
        alleles = String.(group_df[!, allele_sym])
        counts = group_df[!, :count]
        total_count = sum(counts)
        
        # Use existing ratio if available, otherwise calculate
        if "ratio" in names(group_df)
            ratios = group_df[!, :ratio]
            # Sort by ratio (descending) to get major/minor alleles
            sorted_indices = sortperm(ratios, rev=true)
        else
            # Fallback to sorting by count
            sorted_indices = sortperm(counts, rev=true)
        end
        
        sorted_alleles = alleles[sorted_indices]
        sorted_counts = counts[sorted_indices]
        
        # Classify genotype based on number of alleles and ratio
        if length(sorted_alleles) == 1
            genotype_type = "homozygous"
            selected_alleles = [sorted_alleles[1]]
        elseif length(sorted_alleles) == 2
            # Check if minor allele meets ratio threshold
            if "ratio" in names(group_df)
                minor_ratio = ratios[sorted_indices[2]]
                if minor_ratio >= min_ratio
                    genotype_type = "heterozygous"
                    selected_alleles = sorted_alleles[1:2]
                else
                    genotype_type = "homozygous"
                    selected_alleles = [sorted_alleles[1]]
                end
            else
                # Calculate ratio manually
                ratio = sorted_counts[2] / sorted_counts[1]
                if ratio >= min_ratio
                    genotype_type = "heterozygous"
                    selected_alleles = sorted_alleles[1:2]
                else
                    genotype_type = "homozygous"
                    selected_alleles = [sorted_alleles[1]]
                end
            end
        else
            # More than 2 alleles - this is unusual for diploid, mark as uncertain
            genotype_type = "uncertain"
            selected_alleles = sorted_alleles[1:2]  # Take top 2
        end
        
        # Prepare output row
        allele_1 = length(selected_alleles) >= 1 ? selected_alleles[1] : ""
        allele_2 = length(selected_alleles) >= 2 ? selected_alleles[2] : ""
        
        # Get counts for selected alleles (already sorted)
        count_1 = length(sorted_counts) >= 1 ? sorted_counts[1] : 0
        count_2 = length(sorted_counts) >= 2 ? sorted_counts[2] : 0
        
        # Calculate ratio
        ratio = count_1 > 0 ? count_2 / count_1 : 0.0
        
        # Capture other alleles if there are more than 2
        other_alleles = length(sorted_alleles) > 2 ? join(sorted_alleles[3:end], ",") : ""
        
        # Create result tuple with or without novel columns
        if novel_fasta !== nothing
            # Check if alleles are novel
            novel_1 = allele_1 in novel_alleles
            novel_2 = allele_2 in novel_alleles
            
            push!(results, (
                case = case_id,
                gene = gene_id,
                genotype = genotype_type,
                allele_1 = allele_1,
                allele_2 = allele_2,
                count_1 = count_1,
                count_2 = count_2,
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
        uncertain_count = count(x -> x == "uncertain", results_df.genotype)
        
        @info "  Homozygous: $homozygous_count ($(round(100*homozygous_count/total_genotypes, digits=1))%)"
        @info "  Heterozygous: $heterozygous_count ($(round(100*heterozygous_count/total_genotypes, digits=1))%)"
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

end # module
