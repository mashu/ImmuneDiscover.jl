module cooccurrence
    using DataFrames
    using Statistics
    using CSV
    using Distributions
    using Logging
    using Printf
    using Clustering

    export compute_cooccurrence_edges, find_cooccurrence_groups, handle_cooccurrence

    """
        compute_presence_maps(df::DataFrame, case_col::Symbol, allele_col::Symbol)

    Build helpers: unique donors vector, alleles vector, and a Dict allele => Set of donors.
    """
    function compute_presence_maps(df::DataFrame, case_col::Symbol, allele_col::Symbol)
        donors = sort(unique(df[!, case_col]))
        alleles = sort(unique(df[!, allele_col]))
        allele_to_donors = Dict{String, Set{String}}()
        for a in alleles
            # keep donors as a Set for fast intersection
            allele_to_donors[String(a)] = Set{String}(String.(df[df[!, allele_col] .== a, case_col]))
        end
        return donors, alleles, allele_to_donors
    end

    """
        compute_rho_stats(alleles::Vector{T}, allele_to_donors::Dict{String, Set{String}}, N::Int) where {T}

    Compute full pairwise matrices for:
    - rho (phi coefficient)
    - jaccard
    - support (n_shared)
    - p_enrich (hypergeometric one-sided)
    """
    function compute_rho_stats(alleles::Vector, allele_to_donors::Dict{String, Set{String}}, N::Int)
        n = length(alleles)
        R = zeros(Float64, n, n)
        J = zeros(Float64, n, n)
        SUP = zeros(Int, n, n)
        P = ones(Float64, n, n)
        for i in 1:n
            R[i, i] = 1.0
            J[i, i] = 1.0
            SUP[i, i] = 0
            P[i, i] = 0.0
        end
        for i in 1:(n-1)
            ai = String(alleles[i])
            donors_a = allele_to_donors[ai]
            n_a = length(donors_a)
            for j in (i+1):n
                bj = String(alleles[j])
                donors_b = allele_to_donors[bj]
                n_b = length(donors_b)
                # pairwise counts
                n11 = length(intersect(donors_a, donors_b))
                n10 = length(setdiff(donors_a, donors_b))
                n01 = length(setdiff(donors_b, donors_a))
                n00 = N - (n11 + n10 + n01)
                # stats
                ρ = phi_coefficient(n11, n10, n01, n00)
                denom = (n11 + n10 + n01)
                jac = denom == 0 ? 0.0 : n11 / denom
                p = hypergeom_p_enrichment(N, n_a, n_b, n11)
                # assign symmetric
                R[i, j] = ρ; R[j, i] = ρ
                J[i, j] = jac; J[j, i] = jac
                SUP[i, j] = n11; SUP[j, i] = n11
                P[i, j] = p; P[j, i] = p
            end
        end
        return R, J, SUP, P
    end

    """
        print_triangle_diagnostics_matrix(S::Matrix{Float64}, alleles::Vector{String}; threshold::Float64=0.5, max_print::Int=20)

    Triangle diagnostics using similarity matrix S (rho clipped to [0,1]).
    """
    function print_triangle_diagnostics_matrix(S::Matrix{Float64}, alleles::Vector{String}; threshold::Float64=0.5, max_print::Int=20)
        n = length(alleles)
        if n == 0
            @info "Triangle diagnostics: no alleles after filtering"
            return
        end
        # neighbor sets
        nbrs = Dict{Int, Vector{Int}}()
        for i in 1:n
            nbr = Int[]
            for j in 1:n
                if j != i && S[i, j] >= threshold
                    push!(nbr, j)
                end
            end
            nbrs[i] = nbr
        end
        # find triangles
        triangles = Vector{Tuple{Int,Int,Int,Float64,Float64,Float64,Float64}}()
        for i in 1:n
            ni = nbrs[i]
            for a in 1:(length(ni)-1)
                j = ni[a]
                ρij = S[i, j]
                for b in (a+1):length(ni)
                    k = ni[b]
                    ρik = S[i, k]
                    ρjk = S[j, k]
                    if ρjk >= threshold
                        minρ = min(ρij, min(ρik, ρjk))
                        push!(triangles, (i, j, k, ρij, ρik, ρjk, minρ))
                    end
                end
            end
        end
        sort!(triangles, by = x -> x[7], rev=true)
        @info "Triangle diagnostics: alleles=$(n) threshold=$(threshold) triangles=$(length(triangles))"
        for t in triangles[1:min(max_print, length(triangles))]
            @info @sprintf "triangle: (%s, %s, %s) rho=(%.3f, %.3f, %.3f) min_rho=%.3f" alleles[t[1]] alleles[t[2]] alleles[t[3]] t[4] t[5] t[6] t[7]
        end
        if isempty(triangles)
            @info "No triangles at threshold=$(threshold)."
        end
    end

    """
        adjust_bh(pvals::Vector{Float64})

    Benjamini–Hochberg FDR adjustment. Returns q-values aligned to input order.
    """
    function adjust_bh(pvals::Vector{Float64})
        m = length(pvals)
        if m == 0
            return Float64[]
        end
        order = sortperm(pvals)  # indices of p in ascending order
        sorted_p = pvals[order]
        q_sorted = similar(sorted_p)
        min_so_far = 1.0
        for i in reverse(1:m)
            q = (sorted_p[i] * m) / i
            q = q > 1.0 ? 1.0 : q
            if q < min_so_far
                min_so_far = q
            end
            q_sorted[i] = min_so_far
        end
        # map back to original order
        qvals = similar(pvals)
        for i in 1:m
            qvals[order[i]] = q_sorted[i]
        end
        return qvals
    end

    """
        jaccard_index(a::Set{String}, b::Set{String})
    """
    @inline function jaccard_index(a::Set{String}, b::Set{String})
        inter = length(intersect(a, b))
        union_len = length(Base.union(a, b))
        return union_len == 0 ? 0.0 : inter / union_len, inter
    end

    """
        phi_coefficient(n11::Int, n10::Int, n01::Int, n00::Int)

    Phi (rho) coefficient for a 2x2 table.
    """
    @inline function phi_coefficient(n11::Int, n10::Int, n01::Int, n00::Int)
        denom = (n11 + n10) * (n11 + n01) * (n10 + n00) * (n01 + n00)
        if denom <= 0
            return 0.0
        end
        return (n11 * n00 - n10 * n01) / sqrt(denom)
    end

    """
        hypergeom_p_enrichment(N::Int, Ka::Int, Kb::Int, n11::Int)

    One-sided enrichment p-value P[X ≥ n11] where X ~ Hypergeometric(N, Ka, Kb).
    Returns 1.0 if parameters are degenerate.
    """
    function hypergeom_p_enrichment(N::Int, Ka::Int, Kb::Int, n11::Int)
        # Degenerate or independent cases
        if N <= 0 || Ka <= 0 || Kb <= 0 || n11 <= 0
            return 1.0
        end
        # If impossible bounds, return 0 or 1 accordingly
        max_possible = min(Ka, Kb)
        if n11 > max_possible
            return 0.0
        end
        d = Hypergeometric(N, Ka, Kb)
        # survival function: P[X ≥ n11]
        # ccdf is defined as 1 - cdf; ccdf(d, k-1) equals P[X ≥ k]
        return ccdf(d, n11 - 1)
    end

    """
        compute_cooccurrence_edges(df::DataFrame; case_col::AbstractString="case",
                                   allele_col::AbstractString="db_name",
                                   min_donors::Int=1,
                                   min_support::Int=3,
                                   min_jaccard::Float64=0.2)

    Compute pairwise co-occurrence edges using Jaccard index and support (n_shared).
    Returns a DataFrame of qualifying edges and a Dict with per-allele donor sets.
    """
    function compute_cooccurrence_edges(df::DataFrame;
                                        case_col::AbstractString="case",
                                        allele_col::AbstractString="db_name",
                                        min_donors::Int=1,
                                        min_support::Int=3,
                                        min_jaccard::Float64=0.2)
        case_sym = Symbol(case_col)
        allele_sym = Symbol(allele_col)

        # Filter by minimum donors per allele
        allele_counts = combine(groupby(df, allele_sym), nrow => :count)
        valid_alleles = allele_counts[allele_counts.count .>= min_donors, allele_sym]
        filtered_df = filter(x -> x[allele_sym] in valid_alleles, df)

        donors, alleles, allele_to_donors = compute_presence_maps(filtered_df, case_sym, allele_sym)
        N = length(donors)

        rows = Vector{NamedTuple}()
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
                    # 2x2 table components
                    n10 = length(setdiff(donors_a, donors_b))
                    n01 = length(setdiff(donors_b, donors_a))
                    n00 = N - (n11 + n10 + n01)
                    rho = phi_coefficient(n11, n10, n01, n00)
                    p = hypergeom_p_enrichment(N, n_a, n_b, n11)
                    push!(rows, (
                        allele_a = ai,
                        allele_b = bj,
                        n_a = n_a,
                        n_b = n_b,
                        n_shared = n11,
                        jaccard = jac,
                        rho = rho,
                        p_enrich = p
                    ))
                end
            end
        end

        edges = DataFrame(rows)
        # Add BH-FDR q-values
        if nrow(edges) > 0
            qvals = adjust_bh(Vector{Float64}(edges[:, :p_enrich]))
            edges[:, :q_enrich] = qvals
        else
            edges[:, :q_enrich] = Float64[]
        end
        return edges, allele_to_donors
    end

    """
        find_cooccurrence_groups(edges::DataFrame; min_cluster_size::Int=3)

    Build an undirected graph from edges and return connected components (as allele name vectors).
    """
    function find_cooccurrence_groups(edges::DataFrame; min_cluster_size::Int=3)
        if nrow(edges) == 0
            return Vector{Vector{String}}()
        end
        # Build adjacency
        adj = Dict{String, Set{String}}()
        function add_edge(u::String, v::String)
            if u == v
                return
            end
            push!(get!(adj, u, Set{String}()), v)
            push!(get!(adj, v, Set{String}()), u)
        end
        for row in eachrow(edges)
            add_edge(String(row.allele_a), String(row.allele_b))
        end
        # Collect all nodes appearing in edges
        nodes = Set{String}()
        for row in eachrow(edges)
            push!(nodes, String(row.allele_a))
            push!(nodes, String(row.allele_b))
        end
        # Connected components via BFS
        visited = Set{String}()
        components = Vector{Vector{String}}()
        for start in nodes
            if start in visited
                continue
            end
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
            if length(comp) >= min_cluster_size
                push!(components, sort(comp))
            end
        end
        return components
    end

    """
        cluster_hierarchical_from_matrix(S::Matrix{Float64}, alleles::Vector{String};
                                         min_cluster_size::Int=3, threshold::Float64=0.5,
                                         linkage::Symbol=:complete)

    Cluster using hierarchical clustering on a full similarity matrix S (S in [0,1]).
    """
    function cluster_hierarchical_from_matrix(S::Matrix{Float64}, alleles::Vector{String};
                                              min_cluster_size::Int=3, threshold::Float64=0.5,
                                              linkage::Symbol=:complete)
        n = length(alleles)
        if n == 0
            return Vector{Vector{String}}()
        end
        D = 1.0 .- S
        hc = Clustering.hclust(Float64.(D), linkage=linkage)
        labels = Clustering.cutree(hc, h=1.0 - threshold)
        groups = Dict{Int, Vector{String}}()
        for (i, lab) in enumerate(labels)
            push!(get!(groups, lab, String[]), alleles[i])
        end
        comps = Vector{Vector{String}}()
        for g in values(groups)
            if length(g) >= min_cluster_size
                push!(comps, sort(g))
            end
        end
        return comps
    end

    """
        components_from_matrix(S::Matrix{Float64}, alleles::Vector{String}; threshold::Float64=0.5, min_cluster_size::Int=3)

    Build graph by thresholding S and return connected components.
    """
    function components_from_matrix(S::Matrix{Float64}, alleles::Vector{String}; threshold::Float64=0.5, min_cluster_size::Int=3)
        n = length(alleles)
        if n == 0
            return Vector{Vector{String}}()
        end
        adj = Dict{String, Set{String}}()
        for i in 1:n
            for j in (i+1):n
                if S[i, j] >= threshold
                    u = alleles[i]; v = alleles[j]
                    push!(get!(adj, u, Set{String}()), v)
                    push!(get!(adj, v, Set{String}()), u)
                end
            end
        end
        visited = Set{String}()
        components = Vector{Vector{String}}()
        for start in alleles
            if (start in visited) || !(haskey(adj, start))
                continue
            end
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
            if length(comp) >= min_cluster_size
                push!(components, sort(comp))
            end
        end
        return components
    end

    """
        cluster_hierarchical(edges::DataFrame; min_cluster_size::Int=3, threshold::Float64=0.5, linkage::Symbol=:complete)

    Build clusters using hierarchical clustering on rho similarity among alleles present in edges.
    Only pairs present in `edges` contribute similarity; others default to 0.
    Linkage can be :complete, :average, or :single.
    """
    function cluster_hierarchical(edges::DataFrame; min_cluster_size::Int=3, threshold::Float64=0.5, linkage::Symbol=:complete)
        if nrow(edges) == 0
            return Vector{Vector{String}}()
        end
        alleles = sort!(unique(vcat(String.(edges[:, :allele_a]), String.(edges[:, :allele_b]))))
        idx = Dict{String, Int}(a => i for (i, a) in enumerate(alleles))
        n = length(alleles)
        S = zeros(Float64, n, n)
        for i in 1:n
            S[i, i] = 1.0
        end
        for row in eachrow(edges)
            i = idx[String(row.allele_a)]
            j = idx[String(row.allele_b)]
            rho_val = hasproperty(row, :rho) ? row.rho : 0.0
            if rho_val > S[i, j]
                S[i, j] = rho_val
                S[j, i] = rho_val
            end
        end
        D = 1.0 .- S
        hc = Clustering.hclust(Float64.(D), linkage=linkage)
        labels = Clustering.cutree(hc, h=1.0 - threshold)
        groups = Dict{Int, Vector{String}}()
        for (i, lab) in enumerate(labels)
            push!(get!(groups, lab, String[]), alleles[i])
        end
        comps = Vector{Vector{String}}()
        for g in values(groups)
            if length(g) >= min_cluster_size
                push!(comps, sort(g))
            end
        end
        return comps
    end

    """
        handle_cooccurrence(parsed_args, always_gz)

    Handle 'analyze association' command using co-occurrence (Jaccard + support).
    """
    function handle_cooccurrence(parsed_args, always_gz)
        @info "Co-occurrence analysis (Jaccard + support)"
        # Canonical block
        block = parsed_args["analyze"]["cooccurrence"]
        tsv = block["input"]
        case_col = block["case-col"]
        allele_col = block["allele-col"]
        min_donors = block["min-donors"]
        # Similarity/threshold flags are legacy; not used in simplified co-occurrence flow
        min_cluster_size = block["min-cluster-size"]
        cluster_method = get(block, "cluster-method", "components")
        cluster_threshold = get(block, "cluster-threshold", 0.5)
        debug_triangles = get(block, "debug-triangles", false)
        clusters_path = get(block, "clusters", nothing)

        df = CSV.File(tsv, delim='\t') |> DataFrame
        @info "Loaded $(nrow(df)) rows from input file"

        # Filter by min_donors and compute donor presence maps
        case_sym = Symbol(case_col)
        allele_sym = Symbol(allele_col)
        allele_counts = combine(groupby(df, allele_sym), nrow => :count)
        valid_alleles = allele_counts[allele_counts.count .>= min_donors, allele_sym]
        filtered_df = filter(x -> x[allele_sym] in valid_alleles, df)
        donors, alleles, allele_to_donors = compute_presence_maps(filtered_df, case_sym, allele_sym)
        N = length(donors)

        # Compute full rho matrix (plus secondary stats if needed)
        R, J, SUP, P = compute_rho_stats(alleles, allele_to_donors, N)

        if debug_triangles
            S = max.(R, 0.0)
            print_triangle_diagnostics_matrix(S, String.(alleles); threshold=cluster_threshold, max_print=20)
        end

        if clusters_path !== nothing
            # Build similarity matrix from full rho (clip negatives to 0 for similarity)
            S = max.(R, 0.0)
            if cluster_method == "components"
                comps = components_from_matrix(S, String.(alleles); threshold=cluster_threshold, min_cluster_size=min_cluster_size)
            else
                cm = cluster_method == "singel" ? "single" : cluster_method
                link = cm == "average" ? :average : (cm == "single" ? :single : :complete)
                comps = cluster_hierarchical_from_matrix(S, String.(alleles); min_cluster_size=min_cluster_size, threshold=cluster_threshold, linkage=link)
            end
            # Detailed cluster table: one row per allele with its donors
            all_rows = Vector{NamedTuple}()
            for (gid, comp) in enumerate(comps)
                for allele in comp
                    donors_list = sort(collect(get(allele_to_donors, allele, Set{String}())))
                    push!(all_rows, (
                        group_id = gid,
                        allele = allele,
                        donors = join(donors_list, ","),
                        n_donors = length(donors_list)
                    ))
                end
            end
            CSV.write(clusters_path, DataFrame(all_rows), delim='\t')
            @info "Clusters saved in $(clusters_path) ($(length(comps)) groups, $(length(all_rows)) alleles)"
        end
    end

end


