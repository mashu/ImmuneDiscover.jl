module hamming
    using BioAlignments
    using Folds
    using ProgressMeter
    using StringDistances
    using DataFrames
    using MD5
    export hamming_search

    """
        sequence_hash(seq; digits=4)

    Helper function for hashing the sequence
    Note this is not guaranteed to be unique but is simply visual indicator and this encoding is inherited from IgDiscover
    """
    function sequence_hash(seq; digits=4)
        "S"*lpad(string(parse(Int,bytes2hex(MD5.md5(seq))[(end-(digits-1)):end],base=16) % 10^digits), digits, '0')
    end

    """
        unique_name(name, sequence; digits=4)

    Wrapper funciton to update allele names
    """
    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name,"_S")))_$(sequence_hash(sequence,digits=digits))"
    end

    """
        hamming_search(table, db; max_dist=2, column=:genomic_sequence, check_bounds=true, umi=false)

    Function to search demultiplexed reads for sequences that match query with maximum hamming distance
    """
    function hamming_search(table::DataFrame, db::Vector{Tuple{String, String}}; max_dist::Int=2, column::Symbol=:genomic_sequence, check_bounds::Bool=true, umi=false)
        found_list = Vector{Tuple{String, Int, Int, Int, SubString{String}, SubString{String}, SubString{String}, SubString{String}, String, SubString{String}}}()
        @showprogress for subtable in groupby(table, :case)
            case_str = first(subtable.case)
            subtable = subtable[.!occursin.("N", subtable[!, column]), :]
            found = Folds.map(collect(eachrow(subtable))) do row
                process_row(row, db, max_dist, column, case_str, check_bounds, umi=umi)
            end
            append!(found_list, [r for r in found if r !== nothing])
        end
        return found_list
    end
    
    """
        reumi(x)
    
    Helper function to extract UMI from the sequence
    """
    function reumi(x)
        pattern = r"([ACGT]{3}GT[ACGT]{3})([ACGT]{22})([ACGT]{9})$"
        m = match(pattern, x)
        if m === nothing
            return "", "", ""
        end
        return m.captures[1], m.captures[2], m.captures[3]
    end

    """
        process_row(row, db, max_dist, column, case_str, check_bounds; umi=false)

    Helper function to process single row
    """
    function process_row(row::DataFrameRow, db::Vector{Tuple{String, String}}, max_dist::Int, column::Symbol, case_str::String, check_bounds::Bool; umi=false)
        ref = row[column]
        ref_length = length(ref)
        pos = Vector{Tuple{String, Int, Int, Int, SubString{String}}}()
        
        @inbounds for (name, query) in db
            query_length = length(query)
            max_idx = ref_length - query_length + 1
            @inbounds for k in 1:max_idx
                dist = evaluate(Hamming(), query, SubString(ref, k, k + query_length - 1))
                if dist <= max_dist
                    push!(pos, (name, dist, k, k + query_length - 1, SubString(ref, k, k + query_length - 1)))
                end
            end
        end

        if isempty(pos)
            return nothing
        end
        barcode = ""
        if umi
            barcode = reumi(ref)[1]
        end
        best_name, best_dist, start, stop, query_view = first(sort(pos, by=x->x[2]))
        if check_bounds && (ref_length > stop+7) && (start-7 > 0)
            return (best_name, best_dist, start, stop, SubString(ref, start-7, start-1), SubString(ref, start, stop), SubString(ref, stop+1, stop+7), query_view, case_str, barcode)
        elseif !check_bounds
            return (best_name, best_dist, start, stop, "", SubString(ref, start, stop), "", query_view, case_str, barcode)
        end
        return nothing
    end

    """
        summarize(found_list, db; cluster_ratio=0.25, min_count=10)

    Function to summarize results of hamming search
    """
    function summarize(found_list, db; cluster_ratio=0.25, min_count=10)
        lookup = Dict([(y,x) for (x,y) in db])
        found_df = DataFrame(vcat(found_list...))
        rename!(found_df,[:closest_name, :distance,:start,:stop,:prefix,:middle,:suffix,:db_sequence,:case,:barcode])
        sort!(found_df, :closest_name)
        found_df[:,:gene] = [n[1] for n in split.(found_df.closest_name,'*')]
        result = []
        for group in groupby(found_df,[:gene,:case])
            case = string(first(group.case))
            clusters = combine(groupby(group, [:gene, :prefix, :middle, :suffix])) do sdf
                DataFrame(
                    count = nrow(sdf),
                    unique_umi = length(unique(sdf.barcode)),
                    closest_name = unique(sdf.closest_name),
                )
            end
            clusters_filtered = clusters[clusters.count .> (maximum(clusters.count) * cluster_ratio),:]
            clusters_filtered[:,:allele_name] = [get(lookup,seq,unique_name(cname,seq)*" Novel") for (seq, cname) in zip(clusters_filtered.middle, clusters_filtered.closest_name)]
            clusters_filtered[:,:case] .= case
            push!(result, clusters_filtered)
        end
        result_collapsed = vcat(result...)
        result_collapsed_filtered = result_collapsed[result_collapsed.count .>= min_count,:]
        sort!(result_collapsed_filtered, :case)
        result_collapsed_filtered
    end
end
