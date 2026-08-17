# Case / cohort identity for co-occurrence. Dispatch replaces `Union{Dict,Nothing}` maps.

const DEFAULT_POPULATION = "all"

abstract type CaseLabel end
struct BareCase <: CaseLabel
    id::String
end
struct PopCase <: CaseLabel
    id::String
    population::String
end

abstract type Cohort end
struct PooledCohort <: Cohort end
struct StratifiedCohort <: Cohort
    case_to_pop::Dict{String,String}
end

function classify_cases(cases::AbstractVector)
    labeled = PopCase[]
    unlabeled = BareCase[]
    seen = Set{String}()
    for raw in cases
        id = String(raw)
        id in seen && continue
        push!(seen, id)
        parts = split(strip(id), '_')
        if length(parts) < 2
            push!(unlabeled, BareCase(id))
        else
            push!(labeled, PopCase(id, String(parts[end])))
        end
    end
    sort!(labeled; by=c -> c.id)
    sort!(unlabeled; by=c -> c.id)
    return labeled, unlabeled
end

function warn_unlabeled(dropped::AbstractVector{BareCase})
    isempty(dropped) && return
    preview = join((c.id for c in dropped[1:min(end, 5)]), ", ")
    suffix = length(dropped) > 5 ? " ($(length(dropped)) total, showing first 5)" : ""
    @warn "Dropping cases without population suffix in case id: $preview$suffix"
end

function stratified_cohort(cases::AbstractVector)
    labeled, unlabeled = classify_cases(cases)
    isempty(labeled) && error(
        "No cases with population suffix found; name cases as DONOR_POP (e.g. KI_10_EUR)")
    warn_unlabeled(unlabeled)
    return StratifiedCohort(Dict(c.id => c.population for c in labeled)), labeled, unlabeled
end

donor_population(::AbstractString, ::PooledCohort) = DEFAULT_POPULATION
donor_population(case_id::AbstractString, c::StratifiedCohort) = c.case_to_pop[String(case_id)]

function population_breakdown(donors::AbstractVector{<:AbstractString}, ::PooledCohort)
    return DEFAULT_POPULATION, "$(DEFAULT_POPULATION)=$(length(donors))"
end

function population_breakdown(donors::AbstractVector{<:AbstractString}, c::StratifiedCohort)
    counts = Dict{String,Int}()
    for donor in donors
        pop = donor_population(donor, c)
        counts[pop] = get(counts, pop, 0) + 1
    end
    pops = sort(collect(keys(counts)))
    populations = join(pops, ",")
    population_counts = join(["$pop=$(counts[pop])" for pop in pops], ",")
    return populations, population_counts
end

function filter_df_to_cases(df::DataFrame, case_col::Symbol, cases::AbstractSet{String})
    return filter(row -> String(getproperty(row, case_col)) in cases, df)
end

function filter_df_by_population(df::DataFrame, case_col::Symbol,
                                 c::StratifiedCohort, population::AbstractString)
    cases_in_pop = Set(case_id for (case_id, pop) in c.case_to_pop if pop == population)
    return filter_df_to_cases(df, case_col, cases_in_pop)
end

function populations(c::StratifiedCohort)
    return sort(unique(values(c.case_to_pop)))
end
