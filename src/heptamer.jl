module heptamer
    using BioAlignments
    using Folds
    using ProgressMeter
    using StringDistances
    using DataFrames
    using JSON
    include("data.jl")
    using .data

    export extract_heptamers, summarize, load_heptamers

    """
        find_heptamer(suffix, heptamer_iter; max_dist=1, heptamer_length=7)

    Find first occurence of heptamer that is within the hamming distance but prioritize lower distance matches
    """
    function find_heptamer(suffix, heptamer_iter; max_dist=1, heptamer_length=7)
        for d in 0:max_dist # Prioritize smaller distance
            for k in 1:length(suffix)-heptamer_length
                for heptamer in heptamer_iter # Iterate all heptamers in their order
                    if evaluate(Hamming(), heptamer, suffix[k:k+heptamer_length-1]) <= d
                        return k
                    end
                end
            end
        end
        return -1
    end

    """
        extract_heptamers(table, db, heptamers; max_dist=1, b=1, e=8)

    Agument input table with located heptamers and extend trimmed sequence up to that heptamer
    """
    function extract_heptamers(table, db, heptamers; max_dist=1, b=1, e=8)
	@assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
	p = Progress(length(db))
        completed = Folds.map(db) do pair
            query_name, query_seq = pair
            ProgressMeter.next!(p; showvalues = [(:query_name, query_name)])
            query = query_seq[b:end-e]
            subtable = table[occursin.(query, table.genomic_sequence),:]
            short = []
            @inbounds for row in eachrow(subtable)
                genomic_sequence = row[:genomic_sequence]
                interval = findfirst(query, genomic_sequence)
                if interval !== nothing
                    suffix = genomic_sequence[maximum(interval)+1:end]
                    loc = find_heptamer(suffix, heptamers, max_dist=max_dist)
                    if loc > 0
                        extra = suffix[1:loc-1]
                        heptamer = suffix[loc:loc+7-1]
                    else
                        extra = ""
                        heptamer = ""
                    end
                    full_match = ""
                    full_interval = findfirst(query_seq, genomic_sequence)
                    if full_interval !== nothing
                        full_match = string(full_interval)
                    end
                    trimmed_match = string(interval)
                    push!(short, (query_name, full_match, trimmed_match, length(query_seq), length(query*extra), query*extra, heptamer))
                end
            end
            short_df = DataFrame(short)
            if nrow(short_df) > 0
                rename!(short_df, [:db_name, :full_match, :trimmed_match, :db_length, :full_length, :full_sequence, :heptamer])
                sub_short_df = hcat(subtable, short_df)
                sub_short_df
            end
        end
        return vcat([d for d in completed if d !== nothing]...)
    end

    """
        load_heptamers(file)

    Load heptamers from JSON file or create one if not provided.
    """
    function load_heptamers(file)
        if !isfile(file)
            @info "No heptamer file provided, creating one"
            heptamers = Dict("IGKV" => ["CACAGTG","CACACTG","CACTGTG","CACGGTG","CACAATG","CACATTG"],
                             "IGLV" => ["CACAGTG","CACGGTG","CATGGTG","CACGCTG","CACAGCG","CACAGTA","CATAGTG","CACAATG"],
                             "IGHV" => ["CACAGTG","CACAATG","CACAGAG","CACGGTG","CACAGCG"])
            open(file, "w") do f
                JSON.print(f, heptamers)
            end
            @info "Created heptamer JSON $file file"
        end
        if isfile(file)
            heptamers = JSON.parsefile(file)
            @assert all([k ∈ ["IGKV","IGLV","IGHV"] for k in keys(heptamers)]) "JSON file with heptamers must contain IGKV, IGLV and IGHV keys"
        end
        @info "Loaded heptamers for $(join([k for k in keys(heptamers)],','))"
        return heptamers
    end

    """
        summarize(table, query; ratio=0.25, count=10)

    Summarize data by collapsing identical sequences, heptamers and cases
    """
    function summarize(table, query; ratio=0.25, count=10)
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence","db_name","full_match","trimmed_match","db_length","full_length","full_sequence","heptamer"]]) "File must contain following columns well,case,name,genomic_sequence,db_name,full_match,db_length,full_length,full_sequence,heptamer"
        lookup_rev = Dict([(y,x) for (x,y) in query])
        lookup = Dict([(x,y) for (x,y) in query])
        table[:,:trimmed_sequence] = map(x->x.genomic_sequence[UnitRange(parse.(Int,split(x.trimmed_match,':'))...)], eachrow(table))
        table[:,:suffix_sequence] = map(x->x.genomic_sequence[maximum(parse.(Int,split(x.trimmed_match,':'))):end], eachrow(table))
        collapse = []
        for group in groupby(table, [:db_name,:trimmed_sequence,:case])
            full = nrow(group)-sum(ismissing.(group.full_match))
            short = nrow(group)
            for subgroup in groupby(group, [:full_sequence,:heptamer])
                heptamer = first(subgroup.heptamer)
                current_name = string(first(subgroup.db_name))
                new_name = get(lookup_rev, first(subgroup.full_sequence), unique_name(current_name,first(subgroup.full_sequence),digits=4))
                allele_count = nrow(subgroup)
                cases = join(unique(subgroup.case),',')
                full_short_ratio = round(allele_count/short, digits=4)
                full_sequence_length = length(first(subgroup.full_sequence))
                db_sequence_length = length(get(lookup,first(subgroup.db_name),""))
                suffix = first(subgroup.suffix_sequence)
                push!(collapse,
                     (current_name, new_name, full, short, cases, allele_count, full_short_ratio, full_sequence_length, db_sequence_length, first(subgroup.full_sequence), heptamer, suffix))
            end
        end
        collapse_df = DataFrame(collapse)
        rename!(collapse_df, [:db_name, :new_name, :full_count, :trimmed_count, :case, :allele_count, :ratio, :full_length, :db_length, :sequence, :heptamer, :suffix])
        collapse_df_filtered = collapse_df[((collapse_df.ratio .>= ratio) .& (collapse_df.allele_count .>= count)),:]
        unique_df = sort(unique(sort(collapse_df_filtered,:full_count,rev=true), [:new_name, :sequence, :case,:heptamer,:allele_count]),[:case,:new_name])
        unique_df
    end

end
