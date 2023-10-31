module exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds

    """
        exact_search(table, query)

    Exact search per case
    """
    function exact_search(table, query; mincount=10, minfreq=0.2)
        @info "Using mincount: $mincount and minfreq: $minfreq"
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
                    # Comment, because we cannot guarantee thet sequences in database are unique
                    # break
                end
            end
            matches
        end
        result_df = DataFrame(reduce(vcat,[r for r in result if r != []]))
        rename!(result_df, [:case, :db_name])
        final = combine(groupby(result_df,[:case,:db_name]),nrow => :count)
        final[:,:gene] = map(x->first(split(x.db_name,'*')),eachrow(final))
        gdf = groupby(final, [:case, :gene])
        final = transform(gdf, :count => (x->x./maximum(x))=> :frequency)
        filter(r->(r.count >= mincount) & (r.frequency >= minfreq) , final)
    end
end