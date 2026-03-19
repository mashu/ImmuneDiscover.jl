module Cooccurrence
    using DataFrames
    using Statistics
    using CSV
    using Distributions
    using Logging
    using Printf
    using Clustering

    export compute_cooccurrence_edges, find_cooccurrence_groups, handle_cooccurrence
    # Backward compatibility for tests importing `association`
    export compute_edges_and_clusters

    function compute_presence_maps(df::DataFrame, case_col::Symbol, allele_col::Symbol)
        donors = sort(unique(df[!, case_col]))
        alleles = sort(unique(df[!, allele_col]))
        allele_to_donors = Dict{String, Set{String}}()
        for a in alleles
            allele_to_donors[String(a)] = Set{String}(String.(df[df[!, allele_col] .== a, case_col]))
        end
        return donors, alleles, allele_to_donors
    end

    function compute_rho_stats(alleles::Vector, allele_to_donors::Dict{String, Set{String}}, N::Int)
        n = length(alleles)
        R = zeros(Float64, n, n)
        J = zeros(Float64, n, n)
        SUP = zeros(Int, n, n)
        P = ones(Float64, n, n)
        for i in 1:n
            R[i, i] = 1.0; J[i, i] = 1.0; SUP[i, i] = 0; P[i, i] = 0.0
        end
        for i in 1:(n-1)
            ai = String(alleles[i])
            donors_a = allele_to_donors[ai]
            n_a = length(donors_a)
            for j in (i+1):n
                bj = String(alleles[j])
                donors_b = allele_to_donors[bj]
                n_b = length(donors_b)
                n11 = length(intersect(donors_a, donors_b))
                n10 = length(setdiff(donors_a, donors_b))
                n01 = length(setdiff(donors_b, donors_a))
                n00 = N - (n11 + n10 + n01)
                ρ = phi_coefficient(n11, n10, n01, n00)
                denom = (n11 + n10 + n01)
                jac = denom == 0 ? 0.0 : n11 / denom
                p = hypergeom_p_enrichment(N, n_a, n_b, n11)
                R[i, j] = ρ; R[j, i] = ρ
                J[i, j] = jac; J[j, i] = jac
                SUP[i, j] = n11; SUP[j, i] = n11
                P[i, j] = p; P[j, i] = p
            end
        end
        return R, J, SUP, P
    end

    function print_triangle_diagnostics_matrix(S::AbstractMatrix{<:Real}, alleles::AbstractVector{<:AbstractString}; threshold::Float64=0.5, max_print::Int=20)
        n = length(alleles)
        count = 0
        for i in 1:n, j in (i+1):n
            if S[i, j] >= threshold
                count += 1
                if count <= max_print
                    @info @sprintf("  %s — %s: %.3f", alleles[i], alleles[j], S[i, j])
                end
            end
        end
        @info "Total edges above threshold $threshold: $count"
    end

    function adjust_bh(pvals::AbstractVector{<:Real})
        m = length(pvals)
        m == 0 && return Float64[]
        order = sortperm(pvals)
        sorted_p = pvals[order]
        q_sorted = similar(sorted_p)
        min_so_far = 1.0
        for i in reverse(1:m)
            q = min((sorted_p[i] * m) / i, 1.0)
            min_so_far = min(q, min_so_far)
            q_sorted[i] = min_so_far
        end
        qvals = similar(pvals)
        for i in 1:m
            qvals[order[i]] = q_sorted[i]
        end
        return qvals
    end

    @inline function jaccard_index(a::Set{String}, b::Set{String})
        inter = length(intersect(a, b))
        union_len = length(Base.union(a, b))
        return union_len == 0 ? 0.0 : inter / union_len, inter
    end

    @inline function phi_coefficient(n11::Int, n10::Int, n01::Int, n00::Int)
        denom = (n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00)
        denom <= 0 && return 0.0
        return (n11 * n00 - n10 * n01) / sqrt(denom)
    end

    function hypergeom_p_enrichment(N::Int, Ka::Int, Kb::Int, n11::Int)
        (N <= 0 || Ka <= 0 || Kb <= 0 || n11 <= 0) && return 1.0
        n11 > min(Ka, Kb) && return 0.0
        d = Hypergeometric(N, Ka, Kb)
        return ccdf(d, n11 - 1)
    end

    """
        compute_full_stats(df; case_col, allele_col, min_donors)

    Shared computation: filter alleles by donor count, build presence maps,
    compute rho/jaccard/support/pvalue matrices. Returns everything downstream
    code needs (matrices, alleles, allele_to_donors, donor count).
    """
    function compute_full_stats(df::DataFrame;
                                case_col::AbstractString="case",
                                allele_col::AbstractString="db_name",
                                min_donors::Int=1)
        case_sym = Symbol(case_col)
        allele_sym = Symbol(allele_col)
        allele_counts = combine(groupby(df, allele_sym), nrow => :count)
        valid_alleles = allele_counts[allele_counts.count .>= min_donors, allele_sym]
        filtered_df = filter(x -> x[allele_sym] in valid_alleles, df)
        donors, alleles, allele_to_donors = compute_presence_maps(filtered_df, case_sym, allele_sym)
        N = length(donors)
        R, J, SUP, P = compute_rho_stats(alleles, allele_to_donors, N)
        return (; alleles, allele_to_donors, N, R, J, SUP, P)
    end

    function compute_cooccurrence_edges(df::DataFrame;
                                        case_col::AbstractString="case",
                                        allele_col::AbstractString="db_name",
                                        min_donors::Int=1,
                                        min_support::Int=3,
                                        min_jaccard::Float64=0.2)
        stats = compute_full_stats(df; case_col=case_col, allele_col=allele_col, min_donors=min_donors)
        alleles = stats.alleles
        allele_to_donors = stats.allele_to_donors
        N = stats.N
        rows = NamedTuple{(:allele_a,:allele_b,:n_a,:n_b,:n_shared,:jaccard,:rho,:p_enrich), Tuple{String,String,Int,Int,Int,Float64,Float64,Float64}}[]
        for i in 1:(length(alleles)-1)
            ai = String(alleles[i])
            donors_a = allele_to_donors[ai]
            n_a = length(donors_a)
            for j in (i+1):length(alleles)
                bj = String(alleles[j])
                donors_b = allele_to_donors[bj]
                n_b = length(donors_b)
                jac, n11 = jaccard_index(donors_a, donors_b)
                if (n11 >= min_support) && (jac >= min_jaccard)
                    n10 = length(setdiff(donors_a, donors_b))
                    n01 = length(setdiff(donors_b, donors_a))
                    n00 = N - (n11 + n10 + n01)
                    rho = phi_coefficient(n11, n10, n01, n00)
                    p = hypergeom_p_enrichment(N, n_a, n_b, n11)
                    push!(rows, (allele_a=ai, allele_b=bj, n_a=n_a, n_b=n_b, n_shared=n11, jaccard=jac, rho=rho, p_enrich=p))
                end
            end
        end
        edges = DataFrame(rows)
        if nrow(edges) > 0
            edges[:, :q_enrich] = adjust_bh(Vector{Float64}(edges[:, :p_enrich]))
        else
            edges[:, :q_enrich] = Float64[]
        end
        return edges, allele_to_donors
    end

    function find_cooccurrence_groups(edges::DataFrame; min_cluster_size::Int=3)
        nrow(edges) == 0 && return Vector{Vector{String}}()
        adj = Dict{String, Set{String}}()
        for row in eachrow(edges)
            u, v = String(row.allele_a), String(row.allele_b)
            u == v && continue
            push!(get!(adj, u, Set{String}()), v)
            push!(get!(adj, v, Set{String}()), u)
        end
        nodes = Set{String}()
        for row in eachrow(edges)
            push!(nodes, String(row.allele_a))
            push!(nodes, String(row.allele_b))
        end
        visited = Set{String}()
        components = Vector{Vector{String}}()
        for start in nodes
            start in visited && continue
            comp = String[]
            queue = [start]
            push!(visited, start)
            while !isempty(queue)
                u = popfirst!(queue)
                push!(comp, u)
                for v in get(adj, u, Set{String}())
                    if !(v in visited)
                        push!(visited, v)
                        push!(queue, v)
                    end
                end
            end
            length(comp) >= min_cluster_size && push!(components, sort(comp))
        end
        return components
    end

    function cluster_hierarchical_from_matrix(S::AbstractMatrix{<:Real}, alleles::AbstractVector{<:AbstractString};
                                              min_cluster_size::Int=3, threshold::Float64=0.5,
                                              linkage::Symbol=:complete)
        n = length(alleles)
        n == 0 && return Vector{Vector{String}}()
        D = 1.0 .- S
        hc = Clustering.hclust(Float64.(D), linkage=linkage)
        labels = Clustering.cutree(hc, h=1.0 - threshold)
        groups = Dict{Int, Vector{String}}()
        for (i, lab) in enumerate(labels)
            push!(get!(groups, lab, String[]), alleles[i])
        end
        comps = Vector{Vector{String}}()
        for g in values(groups)
            length(g) >= min_cluster_size && push!(comps, sort(g))
        end
        return comps
    end

    function components_from_matrix(S::AbstractMatrix{<:Real}, alleles::AbstractVector{<:AbstractString}; threshold::Float64=0.5, min_cluster_size::Int=3)
        n = length(alleles)
        n == 0 && return Vector{Vector{String}}()
        adj = Dict{String, Set{String}}()
        for i in 1:n, j in (i+1):n
            if S[i, j] >= threshold
                push!(get!(adj, alleles[i], Set{String}()), alleles[j])
                push!(get!(adj, alleles[j], Set{String}()), alleles[i])
            end
        end
        visited = Set{String}()
        components = Vector{Vector{String}}()
        for start in alleles
            (start in visited || !haskey(adj, start)) && continue
            comp = String[]
            queue = [start]
            push!(visited, start)
            while !isempty(queue)
                u = popfirst!(queue)
                push!(comp, u)
                for v in get(adj, u, Set{String}())
                    if !(v in visited)
                        push!(visited, v)
                        push!(queue, v)
                    end
                end
            end
            length(comp) >= min_cluster_size && push!(components, sort(comp))
        end
        return components
    end

    """
        compute_edges_and_clusters(df; kwargs...)

    Backward-compatible wrapper for tests that use `association.compute_edges_and_clusters`.
    """
    function compute_edges_and_clusters(df::DataFrame;
                                        case_col::AbstractString="case",
                                        allele_col::AbstractString="db_name",
                                        min_support::Int=3,
                                        jaccard_threshold::Float64=0.2,
                                        similarity_threshold::Float64=0.5,
                                        min_cluster_size::Int=3)
        edges, allele_to_donors = compute_cooccurrence_edges(df;
            case_col=case_col, allele_col=allele_col,
            min_support=min_support, min_jaccard=jaccard_threshold)
        clusters = find_cooccurrence_groups(edges; min_cluster_size=min_cluster_size)
        # Typed accumulator for detailed cluster rows
        ClusterRow = NamedTuple{(:group_id,:allele,:donors,:n_donors), Tuple{Int,String,String,Int}}
        all_rows = ClusterRow[]
        for (gid, comp) in enumerate(clusters)
            for allele in comp
                donors_list = sort(collect(get(allele_to_donors, allele, Set{String}())))
                push!(all_rows, (group_id=gid, allele=allele, donors=join(donors_list, ","), n_donors=length(donors_list)))
            end
        end
        clusters_detailed = DataFrame(all_rows)
        return edges, clusters, clusters_detailed
    end

    """
        build_edges_from_matrices(R, J, SUP, P, alleles) -> DataFrame

    Convert rho/jaccard/support/pvalue matrices into an edges DataFrame.
    """
    function build_edges_from_matrices(R, J, SUP, P, alleles)
        n = length(alleles)
        EdgeRow = NamedTuple{(:allele_a,:allele_b,:rho,:jaccard,:support,:p_value), Tuple{String,String,Float64,Float64,Int,Float64}}
        rows = EdgeRow[]
        for i in 1:(n-1)
            for j in (i+1):n
                SUP[i,j] > 0 || continue
                push!(rows, (allele_a=String(alleles[i]), allele_b=String(alleles[j]),
                             rho=R[i,j], jaccard=J[i,j], support=SUP[i,j], p_value=P[i,j]))
            end
        end
        return DataFrame(rows)
    end

    function handle_cooccurrence(parsed_args, always_gz)
        @info "Co-occurrence analysis (Jaccard + support)"
        block = parsed_args["analyze"]["cooccurrence"]
        tsv = block["input"]
        case_col = block["case-col"]
        allele_col = block["allele-col"]
        min_donors = block["min-donors"]
        min_cluster_size = block["min-cluster-size"]
        cluster_method = get(block, "cluster-method", "components")
        cluster_threshold = get(block, "cluster-threshold", 0.5)
        debug_triangles = get(block, "debug-triangles", false)
        clusters_path = get(block, "clusters", nothing)

        df = CSV.File(tsv, delim='\t') |> DataFrame
        @info "Loaded $(nrow(df)) rows from input file"

        stats = compute_full_stats(df; case_col=case_col, allele_col=allele_col, min_donors=min_donors)

        edges_df = build_edges_from_matrices(stats.R, stats.J, stats.SUP, stats.P, String.(stats.alleles))
        edges_output = replace(tsv, r"\.(tsv|tsv\.gz)$" => "_edges.tsv")
        if edges_output == tsv
            edges_output = tsv * "_edges.tsv"
        end
        CSV.write(edges_output, edges_df, delim='\t')
        @info "Edges saved to $edges_output ($(nrow(edges_df)) edges from $(length(stats.alleles)) alleles)"

        if debug_triangles
            S = max.(stats.R, 0.0)
            print_triangle_diagnostics_matrix(S, String.(stats.alleles); threshold=cluster_threshold, max_print=20)
        end

        if clusters_path !== nothing
            S = max.(stats.R, 0.0)
            if cluster_method == "components"
                comps = components_from_matrix(S, String.(stats.alleles); threshold=cluster_threshold, min_cluster_size=min_cluster_size)
            else
                link = cluster_method == "average" ? :average : (cluster_method == "single" ? :single : :complete)
                comps = cluster_hierarchical_from_matrix(S, String.(stats.alleles); min_cluster_size=min_cluster_size, threshold=cluster_threshold, linkage=link)
            end
            ClusterRow = NamedTuple{(:group_id,:allele,:donors,:n_donors), Tuple{Int,String,String,Int}}
            all_rows = ClusterRow[]
            for (gid, comp) in enumerate(comps)
                for allele in comp
                    donors_list = sort(collect(get(stats.allele_to_donors, allele, Set{String}())))
                    push!(all_rows, (group_id=gid, allele=allele, donors=join(donors_list, ","), n_donors=length(donors_list)))
                end
            end
            CSV.write(clusters_path, DataFrame(all_rows), delim='\t')
            @info "Clusters saved in $(clusters_path) ($(length(comps)) groups, $(length(all_rows)) alleles)"
        end
    end
end
