module Cooccurrence
    using DataFrames
    using Statistics
    using CSV
    using Distributions
    using Logging
    using Printf
    using Clustering
    using ..Option: Absent, Present, absent, optional, or_default

    export compute_cooccurrence_edges, find_cooccurrence_groups, handle_cooccurrence
    export classify_cases, BareCase, PopCase, PooledCohort, StratifiedCohort
    export build_clusters_dataframe, CompleteLink, AverageLink, SingleLink, ComponentGraph
    export cluster_method

    include("cohort.jl")
    include("cluster_method.jl")
    include("cooccurrence_stats.jl")
    include("cooccurrence_cluster.jl")

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

    function cooccurrence_sidecar_path(tsv::AbstractString, suffix::AbstractString)
        out = replace(tsv, r"\.(tsv|tsv\.gz)$" => suffix)
        return out == tsv ? string(tsv, suffix) : out
    end

    cooccurrence_population_sidecar_path(tsv::AbstractString, population::AbstractString) =
        cooccurrence_sidecar_path(tsv, "_clusters_$(population).tsv")

    """
        build_clusters_dataframe(comps, allele_to_donors, alleles, cohort; analysis_scope)

    One row per allele. Clustered alleles share `group_id` (1-based); unclustered alleles use 0.
    """
    function build_clusters_dataframe(comps::Vector{Vector{String}},
                                      allele_to_donors::Dict{String, Set{String}},
                                      alleles::AbstractVector{<:AbstractString},
                                      cohort::Cohort;
                                      analysis_scope::AbstractString=DEFAULT_POPULATION)
        ClusterRow = NamedTuple{
            (:group_id,:group_size,:analysis_scope,:populations,:population_counts,:allele,:donors,:n_donors),
            Tuple{Int,Int,String,String,String,String,String,Int}}
        rows = ClusterRow[]
        assigned = Set{String}()
        for (gid, comp) in enumerate(comps)
            gsize = length(comp)
            for allele in comp
                push!(assigned, allele)
                donors_list = sort(collect(get(allele_to_donors, allele, Set{String}())))
                pops, pop_counts = population_breakdown(donors_list, cohort)
                push!(rows, (group_id=gid, group_size=gsize, analysis_scope=String(analysis_scope),
                             populations=pops, population_counts=pop_counts, allele=allele,
                             donors=join(donors_list, ","), n_donors=length(donors_list)))
            end
        end
        for allele in alleles
            a = String(allele)
            a in assigned && continue
            donors_list = sort(collect(get(allele_to_donors, a, Set{String}())))
            pops, pop_counts = population_breakdown(donors_list, cohort)
            push!(rows, (group_id=0, group_size=0, analysis_scope=String(analysis_scope),
                         populations=pops, population_counts=pop_counts, allele=a,
                         donors=join(donors_list, ","), n_donors=length(donors_list)))
        end
        clusters_df = DataFrame(rows)
        nrow(clusters_df) == 0 && return clusters_df
        clusters_df[!, :_group_sort] = ifelse.(clusters_df.group_id .== 0, typemax(Int), clusters_df.group_id)
        sort!(clusters_df, [:_group_sort, :allele])
        select!(clusters_df, Not(:_group_sort))
        return clusters_df
    end

    function write_cooccurrence_scope!(df::DataFrame, tsv::AbstractString,
                                       case_col::AbstractString, allele_col::AbstractString,
                                       min_donors::Int, method::ClusterMethod,
                                       cluster_threshold::Float64, min_cluster_size::Int,
                                       analysis_scope::AbstractString,
                                       cohort::Cohort,
                                       clusters_output::AbstractString)
        stats = compute_full_stats(df; case_col=case_col, allele_col=allele_col, min_donors=min_donors)
        length(stats.alleles) == 0 && @warn "No alleles passed filters for analysis_scope=$analysis_scope"
        comps = compute_cluster_components(stats, method;
            cluster_threshold=cluster_threshold, min_cluster_size=min_cluster_size)
        clusters_df = build_clusters_dataframe(comps, stats.allele_to_donors, String.(stats.alleles), cohort;
            analysis_scope=analysis_scope)
        CSV.write(clusters_output, clusters_df, delim='\t')
        n_blocks = length(comps)
        n_block_alleles = count(>(0), clusters_df.group_id)
        @info "Clusters saved to $clusters_output (scope=$analysis_scope, $n_blocks blocks, $n_block_alleles alleles in blocks, $(nrow(clusters_df)) alleles total)"
        return clusters_df
    end

    function write_edges_and_clusters!(df::DataFrame, tsv::AbstractString, case_col::AbstractString,
                                       allele_col::AbstractString, min_donors::Int, method::ClusterMethod,
                                       cluster_threshold::Float64, min_cluster_size::Int,
                                       debug_triangles::Bool, clusters_output::AbstractString,
                                       cohort::Cohort)
        stats = compute_full_stats(df; case_col=case_col, allele_col=allele_col, min_donors=min_donors)
        edges_df = build_edges_from_matrices(stats.R, stats.J, stats.SUP, stats.P, String.(stats.alleles))
        edges_output = cooccurrence_sidecar_path(tsv, "_edges.tsv")
        CSV.write(edges_output, edges_df, delim='\t')
        @info "Edges saved to $edges_output ($(nrow(edges_df)) edges from $(length(stats.alleles)) alleles, scope=$(DEFAULT_POPULATION))"
        if debug_triangles
            S = max.(stats.R, 0.0)
            print_triangle_diagnostics_matrix(S, String.(stats.alleles); threshold=cluster_threshold, max_print=20)
        end
        write_cooccurrence_scope!(df, tsv, case_col, allele_col, min_donors, method,
            cluster_threshold, min_cluster_size, DEFAULT_POPULATION, cohort, clusters_output)
        return stats
    end

    function run_cooccurrence!(::PooledCohort, df::DataFrame, tsv, case_col, allele_col,
                               min_donors, method, cluster_threshold, min_cluster_size,
                               debug_triangles, clusters_output)
        @info "Single-cohort analysis (all cases treated as $(DEFAULT_POPULATION); use --stratify-population for per-population outputs)"
        write_edges_and_clusters!(df, tsv, case_col, allele_col, min_donors, method,
            cluster_threshold, min_cluster_size, debug_triangles, clusters_output, PooledCohort())
        return
    end

    function run_cooccurrence!(cohort::StratifiedCohort, df::DataFrame, tsv, case_col, allele_col,
                               min_donors, method, cluster_threshold, min_cluster_size,
                               debug_triangles, clusters_output)
        case_sym = Symbol(case_col)
        pops = populations(cohort)
        @info "Population stratification on $(length(cohort.case_to_pop)) cases: $(join(pops, ", "))"
        write_edges_and_clusters!(df, tsv, case_col, allele_col, min_donors, method,
            cluster_threshold, min_cluster_size, debug_triangles, clusters_output, cohort)
        for population in pops
            pop_df = filter_df_by_population(df, case_sym, cohort, population)
            nrow(pop_df) == 0 && continue
            pop_output = cooccurrence_population_sidecar_path(tsv, population)
            write_cooccurrence_scope!(pop_df, tsv, case_col, allele_col, min_donors, method,
                cluster_threshold, min_cluster_size, population, cohort, pop_output)
        end
        return
    end

    function handle_cooccurrence(parsed_args)
        @info "Co-occurrence analysis (Jaccard + support)"
        block = parsed_args["analyze"]["cooccurrence"]
        tsv = block["input"]
        case_col = block["case-col"]
        allele_col = block["allele-col"]
        min_donors = block["min-donors"]
        min_cluster_size = block["min-cluster-size"]
        method = cluster_method(get(block, "cluster-method", "complete"))
        cluster_threshold = get(block, "cluster-threshold", 0.7)
        debug_triangles = get(block, "debug-triangles", false)
        stratify_population = get(block, "stratify-population", false)

        df = CSV.File(tsv, delim='\t') |> DataFrame
        @info "Loaded $(nrow(df)) rows from input file"
        clusters_output = or_default(optional(get(block, "clusters", nothing)),
                                     cooccurrence_sidecar_path(tsv, "_clusters.tsv"))

        if stratify_population
            cohort, labeled, _unlabeled = stratified_cohort(df[!, Symbol(case_col)])
            df = filter_df_to_cases(df, Symbol(case_col), Set(c.id for c in labeled))
            run_cooccurrence!(cohort, df, tsv, case_col, allele_col, min_donors, method,
                cluster_threshold, min_cluster_size, debug_triangles, clusters_output)
        else
            run_cooccurrence!(PooledCohort(), df, tsv, case_col, allele_col, min_donors, method,
                cluster_threshold, min_cluster_size, debug_triangles, clusters_output)
        end
    end
end
