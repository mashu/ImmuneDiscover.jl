module linkage
    using DataFrames
    using Statistics
    using SpecialFunctions
    using Printf
    using Clustering
    using CSV

    export compute_linkage_matrices, build_similarity_matrix, perform_clustering, 
           print_clusters, print_top_links, linkage_analysis, linkage_analysis_filtered,
           compute_edges_and_clusters

    """
        compute_linkage_matrices(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol)
    
    Compute linkage disequilibrium (R), Jaccard (J), and support (SUP) matrices from genotype data.
    
    # Arguments
    - `genotypes`: DataFrame with genotype data
    - `case_col`: Symbol for the case/donor column name
    - `allele_col`: Symbol for the allele/db_name column name
    
    # Returns
    - `alleles`: Vector of unique allele names
    - `donors`: Vector of unique donor/case names  
    - `R`: Linkage disequilibrium matrix (Float64 with missings)
    - `J`: Jaccard co-presence matrix
    - `SUP`: Support matrix (co-presence counts)
    """
    function compute_linkage_matrices(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol)
        alleles = sort(unique(genotypes[!, allele_col]))
        donors = sort(unique(genotypes[!, case_col]))
        
        R = zeros(Union{Float64, Missing}, length(alleles), length(alleles))
        J = zeros(Float64, length(alleles), length(alleles))
        SUP = zeros(Int, length(alleles), length(alleles))
        
        for (i, a_allele) in enumerate(alleles)
            for (j, b_allele) in enumerate(alleles)
                if a_allele == b_allele
                    continue
                end
                
                a_cases = genotypes[genotypes[!, allele_col] .== a_allele, case_col]
                b_cases = genotypes[genotypes[!, allele_col] .== b_allele, case_col]

                n11 = length(intersect(a_cases, b_cases)) # both present
                n10 = length(setdiff(a_cases, b_cases)) # a only
                n01 = length(setdiff(b_cases, a_cases)) # b only
                n00 = length(setdiff(donors, a_cases, b_cases)) # neither

                denominator = sqrt((n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00))
                if denominator == 0 || isnan(denominator)
                    r = missing
                else
                    r = (n11 * n00 - n10 * n01) / denominator
                    if isnan(r)
                        r = missing
                    end
                end
                R[i, j] = r
                R[j, i] = r
                
                # Jaccard co-presence (ignore n00)
                denom = n11 + n10 + n01
                jac = denom == 0 ? 0.0 : n11 / denom
                J[i, j] = jac
                J[j, i] = jac
                SUP[i, j] = n11
                SUP[j, i] = n11
            end
        end
        
        return alleles, donors, R, J, SUP
    end

    """
        build_similarity_matrix(R::Matrix{Union{Float64, Missing}}, J::Matrix{Float64}, 
                               SUP::Matrix{Int}; similarity_mode::Symbol=:r, 
                               min_support::Int=3, jaccard_threshold::Float64=0.2)
    
    Build similarity matrix from linkage disequilibrium data with filtering.
    
    # Arguments
    - `R`: Linkage disequilibrium matrix
    - `J`: Jaccard co-presence matrix  
    - `SUP`: Support matrix (co-presence counts)
    - `similarity_mode`: :r2 (LD-style) or :r (positive correlation only)
    - `min_support`: Minimum co-presence count required
    - `jaccard_threshold`: Minimum Jaccard index required
    
    # Returns
    - `S`: Similarity matrix
    """
    function build_similarity_matrix(R::Matrix{Union{Float64, Missing}}, J::Matrix{Float64}, 
                                   SUP::Matrix{Int}; similarity_mode::Symbol=:r, 
                                   min_support::Int=3, jaccard_threshold::Float64=0.2)
        n = size(R, 1)
        S = zeros(Float64, n, n)
        
        for i in 1:n
            for j in 1:n
                if i == j
                    S[i, j] = 1.0
                    continue
                end
                
                rij = R[i, j]
                if ismissing(rij) || (J[i, j] < jaccard_threshold) || (SUP[i, j] < min_support)
                    S[i, j] = 0.0
                else
                    if similarity_mode == :r2
                        S[i, j] = (rij * rij)
                    else
                        S[i, j] = max(rij, 0.0)
                    end
                end
            end
        end
        
        return S
    end

    """
        perform_clustering(S::Matrix{Float64}; similarity_threshold::Float64=0.6)
    
    Perform hierarchical clustering on similarity matrix.
    
    # Arguments
    - `S`: Similarity matrix
    - `similarity_threshold`: Threshold for cutting the dendrogram
    
    # Returns
    - `clusters`: Vector of vectors containing cluster indices
    - `hc`: Hierarchical clustering object
    """
    function perform_clustering(S::Matrix{Float64}; similarity_threshold::Float64=0.6)
        D = 1.0 .- S
        hc = Clustering.hclust(Float64.(D), linkage=:complete)
        distance_threshold = 1.0 - similarity_threshold
        labels = Clustering.cutree(hc, h=distance_threshold)
        
        cluster_map = Dict{Int, Vector{Int}}()
        for (idx, lab) in enumerate(labels)
            push!(get!(cluster_map, lab, Int[]), idx)
        end
        clusters = collect(values(cluster_map))
        
        return clusters, hc
    end

    """
        print_clusters(clusters::Vector{Vector{Int}}, alleles::Vector, S::Matrix{Float64}; 
                      min_cluster_size::Int=3)
    
    Print clusters with statistics.
    
    # Arguments
    - `clusters`: Vector of cluster indices
    - `alleles`: Vector of allele names
    - `S`: Similarity matrix
    - `min_cluster_size`: Minimum cluster size to display
    """
    function print_clusters(clusters::Vector{Vector{Int}}, alleles::Vector, S::Matrix{Float64}; 
                           min_cluster_size::Int=3)
        ck = 0
        for idxs in clusters
            if length(idxs) < min_cluster_size
                continue
            end
            ck += 1
            names = alleles[idxs]
            
            # compute min similarity inside the cluster
            minS = 1.0
            for a in 1:length(idxs)-1
                for b in (a+1):length(idxs)
                    sa = idxs[a]; sb = idxs[b]
                    if S[sa, sb] < minS
                        minS = S[sa, sb]
                    end
                end
            end
            println("Cluster ", ck, " (", length(names), ", minS=", round(minS, digits=3), "): ", join(names, ", "))
        end
    end

    """
        print_top_links(S::Matrix{Float64}, R::Matrix{Union{Float64, Missing}}, 
                       J::Matrix{Float64}, SUP::Matrix{Int}, alleles::Vector; 
                       top_k::Int=30, similarity_mode::Symbol=:r)
    
    Print top pairwise links with statistics.
    
    # Arguments
    - `S`: Similarity matrix
    - `R`: Linkage disequilibrium matrix
    - `J`: Jaccard co-presence matrix
    - `SUP`: Support matrix
    - `alleles`: Vector of allele names
    - `top_k`: Number of top links to display
    - `similarity_mode`: Similarity mode used (for display)
    """
    function print_top_links(S::Matrix{Float64}, R::Matrix{Union{Float64, Missing}}, 
                           J::Matrix{Float64}, SUP::Matrix{Int}, alleles::Vector; 
                           top_k::Int=30, similarity_mode::Symbol=:r)
        pairs = Vector{Tuple{Int,Int,Float64,Float64,Float64,Int}}()
        for i in 1:length(alleles)-1
            for j in (i+1):length(alleles)
                s = S[i, j]
                r = ismissing(R[i, j]) ? 0.0 : R[i, j]
                if s > 0
                    push!(pairs, (i, j, r, s, J[i, j], SUP[i, j]))
                end
            end
        end
        sort!(pairs, by = x -> x[4], rev = true)
        println("Top links (r, S=", String(similarity_mode), ", Jaccard, support):")
        for k in 1:min(top_k, length(pairs))
            (i, j, r, s, jac, sup) = pairs[k]
            println(alleles[i], " -- ", alleles[j], "\tr=", round(r, digits=3), "\tS=", round(s, digits=3), "\tJ=", round(jac, digits=3), "\tn11=", sup)
        end
    end

    """
        linkage_analysis(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol;
                       similarity_mode::Symbol=:r, min_support::Int=3, 
                       jaccard_threshold::Float64=0.2, similarity_threshold::Float64=0.6,
                       min_cluster_size::Int=3, top_k::Int=30, print_results::Bool=true,
                       case_filter_regex::Union{String, Nothing}=nothing)
    
    Complete linkage analysis pipeline.
    
    # Arguments
    - `genotypes`: DataFrame with genotype data
    - `case_col`: Symbol for the case/donor column name
    - `allele_col`: Symbol for the allele/db_name column name
    - `similarity_mode`: :r2 (LD-style) or :r (positive correlation only)
    - `min_support`: Minimum co-presence count required
    - `jaccard_threshold`: Minimum Jaccard index required
    - `similarity_threshold`: Threshold for cutting the dendrogram
    - `min_cluster_size`: Minimum cluster size to display
    - `top_k`: Number of top links to display
    - `print_results`: Whether to print cluster and link results
    - `case_filter_regex`: Regex pattern to filter cases (default: "[ACDERF]" for monkey cases)
    
    # Returns
    - `alleles`: Vector of unique allele names
    - `donors`: Vector of unique donor/case names
    - `R`: Linkage disequilibrium matrix
    - `J`: Jaccard co-presence matrix
    - `SUP`: Support matrix
    - `S`: Similarity matrix
    - `clusters`: Vector of cluster indices
    - `hc`: Hierarchical clustering object
    """
    function linkage_analysis(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol;
                             similarity_mode::Symbol=:r, min_support::Int=3, 
                             jaccard_threshold::Float64=0.2, similarity_threshold::Float64=0.6,
                             min_cluster_size::Int=3, top_k::Int=30, print_results::Bool=true,
                             case_filter_regex::Union{String, Nothing}=nothing)
        
        # Apply case filtering if regex is provided
        filtered_genotypes = copy(genotypes)
        if case_filter_regex !== nothing
            filter!(x -> startswith(x[case_col], Regex(case_filter_regex)), filtered_genotypes)
        end
        
        # Compute linkage matrices
        alleles, donors, R, J, SUP = compute_linkage_matrices(filtered_genotypes, case_col, allele_col)
        
        # Build similarity matrix
        S = build_similarity_matrix(R, J, SUP; 
                                  similarity_mode=similarity_mode, 
                                  min_support=min_support, 
                                  jaccard_threshold=jaccard_threshold)
        
        # Perform clustering
        clusters, hc = perform_clustering(S; similarity_threshold=similarity_threshold)
        
        # Print results if requested
        if print_results
            print_clusters(clusters, alleles, S; min_cluster_size=min_cluster_size)
            print_top_links(S, R, J, SUP, alleles; top_k=top_k, similarity_mode=similarity_mode)
        end
        
        return alleles, donors, R, J, SUP, S, clusters, hc
    end

    """
        compute_edges_and_clusters(df::DataFrame; case_col::AbstractString="case", 
                                  allele_col::AbstractString="db_name", min_donors::Int=1,
                                  min_support::Int=3, jaccard_threshold::Float64=0.2,
                                  similarity::Symbol=:r, similarity_threshold::Float64=0.5,
                                  min_cluster_size::Int=3)
    
    Compute pairwise metrics table (edges) and complete-linkage clusters using positive-only r similarity 
    gated by co-presence (Jaccard) and support (n11). Returns (edges_df, clusters, clusters_detailed).
    
    This function provides compatibility with the existing CLI interface.
    Note: Case filtering should be applied at the demultiplex stage, not here.
    
    # Arguments
    - `df`: DataFrame with genotype data
    - `case_col`: Name of the case/donor column
    - `allele_col`: Name of the allele column
    - `min_donors`: Minimum donors required to include an allele
    - `min_support`: Minimum co-presence count required
    - `jaccard_threshold`: Minimum Jaccard index required
    - `similarity`: Similarity mode (:r or :r2)
    - `similarity_threshold`: Threshold for cutting the dendrogram
    - `min_cluster_size`: Minimum cluster size to output
    
    # Returns
    - `edges_df`: DataFrame with pairwise metrics
    - `clusters`: Vector of vectors containing allele names in each cluster
    - `clusters_detailed`: Vector of DataFrames with detailed cluster information including:
        - allele: allele name
        - donors: comma-separated list of donors having this allele
        - n_donors: number of donors having this allele
        - mean_n11, mean_n10, mean_n01, mean_n00: average pairwise counts
        - mean_r, max_r, min_r: average, maximum, and minimum linkage disequilibrium r values
    """
    function compute_edges_and_clusters(df::DataFrame; 
                                      case_col::AbstractString="case", 
                                      allele_col::AbstractString="db_name", 
                                      min_donors::Int=1,
                                      min_support::Int=3, 
                                      jaccard_threshold::Float64=0.2,
                                      similarity::Symbol=:r, 
                                      similarity_threshold::Float64=0.5,
                                      min_cluster_size::Int=3)
        
        # Convert string column names to symbols
        case_sym = Symbol(case_col)
        allele_sym = Symbol(allele_col)
        
        # Filter alleles by minimum donor count
        allele_counts = combine(groupby(df, allele_sym), nrow => :count)
        valid_alleles = allele_counts[allele_counts.count .>= min_donors, allele_sym]
        filtered_df = filter(x -> x[allele_sym] in valid_alleles, df)
        
        # Run linkage analysis
        alleles, donors, R, J, SUP, S, clusters, hc = linkage_analysis(
            filtered_df, case_sym, allele_sym;
            similarity_mode=similarity,
            min_support=min_support,
            jaccard_threshold=jaccard_threshold,
            similarity_threshold=similarity_threshold,
            min_cluster_size=min_cluster_size,
            print_results=false
        )
        
        # Create edges DataFrame
        edges_rows = Vector{NamedTuple}()
        for i in 1:length(alleles)-1
            for j in (i+1):length(alleles)
                r_val = ismissing(R[i, j]) ? 0.0 : R[i, j]
                jaccard = J[i, j]
                support = SUP[i, j]
                similarity_val = S[i, j]
                
                if similarity_val > 0  # Only include edges with positive similarity
                    push!(edges_rows, (
                        allele_a = alleles[i],
                        allele_b = alleles[j],
                        r = r_val,
                        r2 = r_val * r_val,
                        jaccard = jaccard,
                        support = support,
                        similarity = similarity_val
                    ))
                end
            end
        end
        
        edges_df = DataFrame(edges_rows)
        
        # Convert cluster indices to allele names and collect detailed information
        clusters_alleles = Vector{Vector{String}}()
        clusters_detailed = Vector{DataFrame}()
        
        for cluster in clusters
            if length(cluster) >= min_cluster_size
                cluster_alleles = String.(alleles[cluster])
                push!(clusters_alleles, cluster_alleles)
                
                # Create detailed cluster information
                cluster_rows = Vector{NamedTuple}()
                for allele_idx in cluster
                    allele_name = String(alleles[allele_idx])
                    
                    # Get donors for this allele
                    allele_donors = filtered_df[filtered_df[!, allele_sym] .== alleles[allele_idx], case_sym]
                    donors_list = sort(unique(allele_donors))
                    donors_str = join(donors_list, ",")
                    
                    # Calculate pairwise statistics with other alleles in cluster
                    n11_values = Float64[]
                    n10_values = Float64[]
                    n01_values = Float64[]
                    n00_values = Float64[]
                    r_values = Float64[]
                    
                    for other_idx in cluster
                        if other_idx != allele_idx
                            other_allele = alleles[other_idx]
                            other_donors = filtered_df[filtered_df[!, allele_sym] .== other_allele, case_sym]
                            
                            n11 = length(intersect(allele_donors, other_donors))
                            n10 = length(setdiff(allele_donors, other_donors))
                            n01 = length(setdiff(other_donors, allele_donors))
                            n00 = length(setdiff(donors, allele_donors, other_donors))
                            
                            r = (n11 * n00 - n10 * n01) / sqrt((n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00))
                            if isnan(r)
                                r = 0.0
                            end
                            
                            push!(n11_values, n11)
                            push!(n10_values, n10)
                            push!(n01_values, n01)
                            push!(n00_values, n00)
                            push!(r_values, r)
                        end
                    end
                    
                    # Calculate summary statistics
                    mean_n11 = length(n11_values) > 0 ? mean(n11_values) : 0.0
                    mean_n10 = length(n10_values) > 0 ? mean(n10_values) : 0.0
                    mean_n01 = length(n01_values) > 0 ? mean(n01_values) : 0.0
                    mean_n00 = length(n00_values) > 0 ? mean(n00_values) : 0.0
                    mean_r = length(r_values) > 0 ? mean(r_values) : 0.0
                    max_r = length(r_values) > 0 ? maximum(r_values) : 0.0
                    min_r = length(r_values) > 0 ? minimum(r_values) : 0.0
                    
                    push!(cluster_rows, (
                        allele = allele_name,
                        donors = donors_str,
                        n_donors = length(donors_list),
                        mean_n11 = round(mean_n11, digits=3),
                        mean_n10 = round(mean_n10, digits=3),
                        mean_n01 = round(mean_n01, digits=3),
                        mean_n00 = round(mean_n00, digits=3),
                        mean_r = round(mean_r, digits=3),
                        max_r = round(max_r, digits=3),
                        min_r = round(min_r, digits=3)
                    ))
                end
                
                cluster_df = DataFrame(cluster_rows)
                push!(clusters_detailed, cluster_df)
            end
        end
        
        return edges_df, clusters_alleles, clusters_detailed
    end

    """
        linkage_analysis_filtered(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol;
                                 similarity_mode::Symbol=:r, min_support::Int=3, 
                                 jaccard_threshold::Float64=0.2, similarity_threshold::Float64=0.6,
                                 min_cluster_size::Int=3, top_k::Int=30, print_results::Bool=true,
                                 case_filter_regex::String="[ACDERF]")
    
    Convenience function for linkage analysis with case filtering.
    Equivalent to calling linkage_analysis with the specified case_filter_regex.
    
    # Arguments
    - `genotypes`: DataFrame with genotype data
    - `case_col`: Symbol for the case/donor column name
    - `allele_col`: Symbol for the allele/db_name column name
    - `similarity_mode`: :r2 (LD-style) or :r (positive correlation only)
    - `min_support`: Minimum co-presence count required
    - `jaccard_threshold`: Minimum Jaccard index required
    - `similarity_threshold`: Threshold for cutting the dendrogram
    - `min_cluster_size`: Minimum cluster size to display
    - `top_k`: Number of top links to display
    - `print_results`: Whether to print cluster and link results
    - `case_filter_regex`: Regex pattern to filter cases (default: "[ACDERF]")
    
    # Returns
    - `alleles`: Vector of unique allele names
    - `donors`: Vector of unique donor/case names
    - `R`: Linkage disequilibrium matrix
    - `J`: Jaccard co-presence matrix
    - `SUP`: Support matrix
    - `S`: Similarity matrix
    - `clusters`: Vector of cluster indices
    - `hc`: Hierarchical clustering object
    """
    function linkage_analysis_filtered(genotypes::DataFrame, case_col::Symbol, allele_col::Symbol;
                                      similarity_mode::Symbol=:r, min_support::Int=3, 
                                      jaccard_threshold::Float64=0.2, similarity_threshold::Float64=0.6,
                                      min_cluster_size::Int=3, top_k::Int=30, print_results::Bool=true,
                                      case_filter_regex::String="[ACDERF]")
        return linkage_analysis(genotypes, case_col, allele_col;
                              similarity_mode=similarity_mode, min_support=min_support,
                              jaccard_threshold=jaccard_threshold, similarity_threshold=similarity_threshold,
                              min_cluster_size=min_cluster_size, top_k=top_k, print_results=print_results,
                              case_filter_regex=case_filter_regex)
    end

end