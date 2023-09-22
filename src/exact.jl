module exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds

    """
        exact_search(table, query)

    Exact search per case
    """
    function exact_search(table, query)
        @assert all(names(table) .== ["case","name","genomic_sequence"]) "File must contain following columns case,name,genomic_sequence"
        p = Progress(nrow(table))
        result = Folds.map(eachrow(table)) do row
            next!(p)
            case = row.case
            matches = Vector{Tuple}()
            @inbounds for (name, seq) in query
                match = occursin(seq, row.genomic_sequence)
                if match
                    push!(matches, (string(case), string(name)))
                    break
                end
            end
            matches
        end
        result_df = DataFrame(reduce(vcat,[r for r in result if r != []]))
        rename!(result_df, [:case, :db_name])
        combine(groupby(result_df,[:case,:db_name]),nrow => :count)
    end
end