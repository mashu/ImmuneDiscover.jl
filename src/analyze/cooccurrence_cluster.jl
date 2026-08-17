function adjacency_from_edges(edges::DataFrame)
    adj = Dict{String, Set{String}}()
    for row in eachrow(edges)
        u, v = String(row.allele_a), String(row.allele_b)
        u == v && continue
        push!(get!(adj, u, Set{String}()), v)
        push!(get!(adj, v, Set{String}()), u)
    end
    return adj
end

function adjacency_from_matrix(S::AbstractMatrix{<:Real}, alleles::AbstractVector{<:AbstractString}, threshold::Float64)
    n = length(alleles)
    adj = Dict{String, Set{String}}()
    for i in 1:n, j in (i+1):n
        S[i, j] >= threshold || continue
        push!(get!(adj, alleles[i], Set{String}()), alleles[j])
        push!(get!(adj, alleles[j], Set{String}()), alleles[i])
    end
    return adj
end

function connected_components(adj::Dict{String,Set{String}}, nodes; min_cluster_size::Int)
    visited = Set{String}()
    components = Vector{Vector{String}}()
    for start in nodes
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

function find_cooccurrence_groups(edges::DataFrame; min_cluster_size::Int=3)
    nrow(edges) == 0 && return Vector{Vector{String}}()
    adj = adjacency_from_edges(edges)
    nodes = Set{String}()
    for row in eachrow(edges)
        push!(nodes, String(row.allele_a))
        push!(nodes, String(row.allele_b))
    end
    return connected_components(adj, nodes; min_cluster_size=min_cluster_size)
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

function components_from_matrix(S::AbstractMatrix{<:Real}, alleles::AbstractVector{<:AbstractString};
                                threshold::Float64=0.5, min_cluster_size::Int=3)
    n = length(alleles)
    n == 0 && return Vector{Vector{String}}()
    return connected_components(adjacency_from_matrix(S, alleles, threshold), alleles;
                                min_cluster_size=min_cluster_size)
end
