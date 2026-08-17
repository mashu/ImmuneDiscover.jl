# Clustering on the co-occurrence similarity matrix. One type per method so the handler
# parses the CLI string once and everything after that is dispatch.

abstract type ClusterMethod end
struct CompleteLink <: ClusterMethod end
struct AverageLink <: ClusterMethod end
struct SingleLink <: ClusterMethod end
struct ComponentGraph <: ClusterMethod end

cluster_method(name::AbstractString) = cluster_method(Val{Symbol(name)}())
cluster_method(::Val{:complete}) = CompleteLink()
cluster_method(::Val{:average}) = AverageLink()
cluster_method(::Val{:single}) = SingleLink()
cluster_method(::Val{:components}) = ComponentGraph()
cluster_method(::Val{S}) where {S} = error("Unknown cluster method: $S")

hclust_linkage(::CompleteLink) = :complete
hclust_linkage(::AverageLink) = :average
hclust_linkage(::SingleLink) = :single

function compute_cluster_components(stats, method::ComponentGraph;
                                    cluster_threshold::Float64=0.7, min_cluster_size::Int=3)
    S = max.(stats.R, 0.0)
    return components_from_matrix(S, String.(stats.alleles);
        threshold=cluster_threshold, min_cluster_size=min_cluster_size)
end

function compute_cluster_components(stats, method::ClusterMethod;
                                    cluster_threshold::Float64=0.7, min_cluster_size::Int=3)
    S = max.(stats.R, 0.0)
    return cluster_hierarchical_from_matrix(S, String.(stats.alleles);
        min_cluster_size=min_cluster_size, threshold=cluster_threshold,
        linkage=hclust_linkage(method))
end
