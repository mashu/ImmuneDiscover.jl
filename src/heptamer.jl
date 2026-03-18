module heptamer
    using BioAlignments
    using Folds
    using ProgressMeter
    using StringDistances
    using DataFrames
    using JSON
    using CSV
    using ..data

    export extract_heptamers, summarize, load_heptamers, handle_heptamer

    function find_heptamer(suffix, heptamer_iter; max_dist=1, heptamer_length=7)
        for d in 0:max_dist
            for k in 1:length(suffix)-heptamer_length
                for heptamer in heptamer_iter
                    if evaluate(Hamming(), heptamer, suffix[k:k+heptamer_length-1]) <= d
                        return k
                    end
                end
            end
        end
        return -1
    end

    function extract_heptamers(table, db, heptamers; max_dist=1, b=1, e=8)
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        p = Progress(length(db))
        completed = Folds.map(db) do pair
            query_name, query_seq = pair
            ProgressMeter.next!(p; showvalues = [(:query_name, query_name)])
            query = query_seq[b:end-e]
            subtable = table[occursin.(query, table.genomic_sequence),:]
            short = Vector{Tuple{String, String, String, Int, Int, String, String}}()
            @inbounds for row in eachrow(subtable)
                genomic_sequence = row[:genomic_sequence]
                interval = findfirst(query, genomic_sequence)
                if interval !== nothing
                    suffix = genomic_sequence[maximum(interval)+1:end]
                    loc = find_heptamer(suffix, heptamers, max_dist=max_dist)
                    if loc > 0
                        extra = suffix[1:loc-1]
                        heptamer_seq = suffix[loc:loc+7-1]
                    else
                        extra = ""
                        heptamer_seq = ""
                    end
                    full_match = ""
                    full_interval = findfirst(query_seq, genomic_sequence)
                    if full_interval !== nothing
                        full_match = string(full_interval)
                    end
                    trimmed_match = string(interval)
                    push!(short, (query_name, full_match, trimmed_match, length(query_seq), length(query*extra), query*extra, heptamer_seq))
                end
            end
            short_df = DataFrame(short)
            if nrow(short_df) > 0
                rename!(short_df, [:db_name, :full_match, :trimmed_match, :db_length, :full_length, :sequence, :heptamer])
            end
            short_df
        end
        result = vcat(completed...)
        return result
    end

    function summarize(df, db; ratio=0.2, count=10)
        nrow(df) == 0 && return DataFrame()
        summary = combine(groupby(df, [:db_name, :sequence, :heptamer]), nrow => :count)
        transform!(groupby(summary, :db_name), :count => (x -> x ./ maximum(x)) => :ratio)
        filter!(x -> x.count >= count && x.ratio >= ratio, summary)
        sort!(summary, [:db_name, :count], rev=[false, true])
        return summary
    end

    function load_heptamers(json_path)
        return JSON.parsefile(json_path)
    end

    function handle_heptamer(parsed_args, immunediscover_module, always_gz)
        @info "Heptamer search"
        table = immunediscover_module.load_demultiplex(parsed_args["search"]["heptamer"]["tsv"])
        db = immunediscover_module.load_fasta(parsed_args["search"]["heptamer"]["fasta"], validate=false)
        chain = parsed_args["search"]["heptamer"]["chain"]
        heptamers = load_heptamers(parsed_args["search"]["heptamer"]["json"])
        heptamer_df = extract_heptamers(table, db, heptamers[chain];
            max_dist=parsed_args["search"]["heptamer"]["maxdist"],
            b=parsed_args["search"]["heptamer"]["begin"]+1,
            e=parsed_args["search"]["heptamer"]["end"])
        output = always_gz(parsed_args["search"]["heptamer"]["output"])
        CSV.write(output, heptamer_df, compress=true, delim='\t')
        @info "Heptamer search results saved to $output"
        summary_df = summarize(heptamer_df, db,
            ratio=parsed_args["search"]["heptamer"]["ratio"],
            count=parsed_args["search"]["heptamer"]["mincount"])
        summary_file = parsed_args["search"]["heptamer"]["summary"]
        CSV.write(summary_file, summary_df, delim='\t')
        @info "Heptamer summary saved to $summary_file"
    end
end
