module hamming
    using BioAlignments
    using Folds
    using ProgressMeter
    using StringDistances
    using DataFrames
    export hamming_search

    """
    Helper function for hashing the sequence
    Note this is not guaranteed to be unique but is simply visual indicator and this encoding is inherited from IgDiscover
    """
    function sequence_hash(seq; digits=4)
        "S"*lpad(string(parse(Int,bytes2hex(MD5.md5(seq))[(end-(digits-1)):end],base=16) % 10^digits), digits, '0')
    end

    """
    Wrapper funciton to update allele names
    """
    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name,"_S")))_$(sequence_hash(sequence,digits=digits))"
    end

    """
    Function to search demultiplexed reads for sequences that match query with maximum hamming distance
    """
    function hamming_search(table, db;max_dist=2)
        found_list = []
        @showprogress for subtable in groupby(table,:case)
            case = string(first(subtable.case))
            subtable = subtable[occursin.("N",subtable.genomic_sequence) .!= 1,:]
            rows = Folds.map(collect(eachrow(subtable))) do row
                case = string(row.case)
                ref = row.genomic_sequence
                pos = []
                ref_length = length(ref)
                for (name, query) in db
                    query_length = length(query)
                    for k in 1:length(ref)-query_length
                        if (k+query_length-1) < ref_length
                            dist = evaluate(Hamming(),query, ref[k:k+query_length-1])
                            if dist <= max_dist
                                push!(pos, (name, dist, k,k+query_length-1, query))
                            end
                        end
                    end
                end
                if length(pos) > 0
                    best_name, best_dist, start, stop, query = first(sort(pos,by=x->x[2]))    
                    if (ref_length > stop+7) & ((start-7) > 0)
                        match_sequence = ref[start:stop]
                        match_prefix = ref[start-7:start-1]
                        match_suffix = ref[stop+1:stop+7]
                        (best_name, best_dist, start, stop, match_prefix, match_sequence, match_suffix, query, case)
                    end
                end
            end
            found = [r for r in rows if r != nothing]
            push!(found_list,found)
        end
        return found_list
    end

    function summarize(found_list, db; cluster_ratio=0.25, min_count=10)
        lookup = Dict([(y,x) for (x,y) in db])
        found_df = DataFrame(vcat(found_list...))
        rename!(found_df,[:closest_name, :distance,:start,:stop,:prefix,:middle,:suffix,:db_sequence,:case])
        sort!(found_df, :closest_name)
        found_df[:,:gene] = [n[1] for n in split.(found_df.closest_name,'*')]
        result = []
        for group in groupby(found_df,[:gene,:case])
            case = string(first(group.case))
            clusters = combine(groupby(group,[:gene,:prefix,:middle,:suffix]),nrow => :cluster_size)
            clusters_filtered = clusters[clusters.cluster_size .> (maximum(clusters.cluster_size) * cluster_ratio),:]
            clusters_filtered[:,:allele_name] = [get(lookup,seq,immune.unique_name(string(first(group.gene)),seq)*" Novel") for seq in clusters_filtered.middle]
            clusters_filtered[:,:case] .= case
            push!(result, clusters_filtered)
        end
        result_collapsed = vcat(result...)
        result_collapsed_filtered = result_collapsed[result_collapsed.cluster_size .>= min_count,:]
        sort!(result_collapsed_filtered, :case)
        result_collapsed_filtered
    end
end