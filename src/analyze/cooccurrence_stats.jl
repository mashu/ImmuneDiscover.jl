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
compute rho/jaccard/support/pvalue matrices.
"""
function compute_full_stats(df::DataFrame;
                            case_col::AbstractString="case",
                            allele_col::AbstractString="db_name",
                            min_donors::Int=1)
    case_sym = Symbol(case_col)
    allele_sym = Symbol(allele_col)
    donor_counts = combine(groupby(df, allele_sym), case_sym => (x -> length(unique(x))) => :n_donors)
    valid_alleles = Set(donor_counts[donor_counts.n_donors .>= min_donors, allele_sym])
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
    R, J, SUP, P = stats.R, stats.J, stats.SUP, stats.P
    rows = NamedTuple{(:allele_a,:allele_b,:n_a,:n_b,:n_shared,:jaccard,:rho,:p_enrich), Tuple{String,String,Int,Int,Int,Float64,Float64,Float64}}[]
    for i in 1:(length(alleles)-1)
        ai = String(alleles[i])
        n_a = length(allele_to_donors[ai])
        for j in (i+1):length(alleles)
            n11 = SUP[i, j]
            jac = J[i, j]
            if (n11 >= min_support) && (jac >= min_jaccard)
                bj = String(alleles[j])
                push!(rows, (allele_a=ai, allele_b=bj, n_a=n_a, n_b=length(allele_to_donors[bj]),
                             n_shared=n11, jaccard=jac, rho=R[i, j], p_enrich=P[i, j]))
            end
        end
    end
    edges = DataFrame(rows)
    edges[:, :q_enrich] = nrow(edges) > 0 ? adjust_bh(Vector{Float64}(edges[:, :p_enrich])) : Float64[]
    return edges, allele_to_donors
end

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
    edges = DataFrame(rows)
    edges[:, :q_value] = nrow(edges) > 0 ? adjust_bh(Vector{Float64}(edges[:, :p_value])) : Float64[]
    return edges
end
