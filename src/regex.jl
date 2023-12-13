module regex
    using ProgressMeter
    using DataFrames
    using Folds
    using StringDistances
    using BioAlignments
    using DataStructures
    include("data.jl")
    using .data
    export search_regex, search_frequent_flanks, annotate, update_names!, flanks, longest_common_substring, combine_accumulators

    """
        longest_common_substring(str1::String, str2::String)::String

    Find longest common substring between two strings
    """
    function longest_common_substring(str1::String, str2::String)::String
        len1, len2 = length(str1), length(str2)
        dp = zeros(Int, len1 + 1, len2 + 1)
        
        longest_len = 0
        end_index = 0
    
        for i in 1:len1
            for j in 1:len2
                if str1[i] == str2[j]
                    dp[i + 1, j + 1] = dp[i, j] + 1
                    if dp[i + 1, j + 1] > longest_len
                        longest_len = dp[i + 1, j + 1]
                        end_index = i
                    end
                else
                    dp[i + 1, j + 1] = 0
                end
            end
        end
    
        if longest_len == 0
            return ""
        else
            return str1[end_index - longest_len + 1:end_index]
        end
    end

    """
    Extract flanks on both sides
    """
    function flanks(query, target; nprefix=7, nsuffix=7)
        range = findfirst(query, target)
        return target[minimum(range)-nprefix:minimum(range)-1],target[maximum(range)+1:maximum(range)+nsuffix]
    end

    """
        update_names!(output_df)
    
    Update names of alleles if they're novel given the distance
    """
    function update_names!(output_df)
        for i in 1:nrow(output_df)
            if output_df[i, :distance] > 0
                output_df[i, :best_name] = unique_name(output_df[i,:best_name], output_df[i,:best_aln]) * " Novel"
            end
        end
    end

    """
        combine_accumulators(accumulators::Vector{Accumulator{String,Int}})
    
    Combine accumulators
    """
    function combine_accumulators(accumulators::Vector{Accumulator{String,Int}})
        # Create a master accumulator
        master_accumulator = Accumulator{String, Int}()
    
        # Iterate through each accumulator
        for acc in accumulators
            # Iterate through each element in the accumulator
            for (key, count) in acc
                # Add the count to the master accumulator
                master_accumulator[key] += count
            end
        end
    
        return master_accumulator
    end

    """
    Collect frequent flanks
    """
    function search_frequent_flanks(table, db; min_length=16, min_frequency=0.5, min_count=100, nprefix=7, nsuffix=7)
        p = Progress(nrow(table), desc = "Collecting prefix/suffix")
        @info "Prefix and suffix frequency threshold set to $min_frequency and count threshold set to $min_count"
        @info "Prefixes are set to be $(nprefix) bases long and suffixes must be $(nsuffix) bases long"
        result = Folds.map(db) do (dname, dseq)
            prefix_set = Accumulator{String,Int}()
            suffix_set = Accumulator{String,Int}()
            if length(dseq) < min_length
                @warn "$dname shorter than 16, skipping"
                return
            end
            subdf = table[occursin.(string(dseq), table.genomic_sequence),:]
            hits = transform(subdf, :genomic_sequence => ByRow(s-> flanks(string(dseq),s, nprefix=nprefix, nsuffix=nsuffix)) => [:prefix,:suffix])
            combined = combine(groupby(hits,[:prefix,:suffix, :well, :case]),nrow => :count)
            if nrow(combined) == 0
                return
            end
            # Prefixes
            frequent = combined[(combined.count ./ maximum(combined.count)) .>= min_frequency,:]
            if maximum(combined.count) > min_count
                for p in frequent.prefix
                    push!(prefix_set, p)
                end
                for s in frequent.suffix
                    push!(suffix_set, s)
                end
            end
            next!(p)
            return (dname=dname, prefix=prefix_set, suffix=suffix_set)
        end
        finish!(p)
        return filter(x->x !== nothing, result)
    end

    """
        create_regex(prefix_sorted, suffix_sorted, minlen, maxlen)

    Create regex from frequent flanks
    """
    function create_regex(prefix_sorted, suffix_sorted, minlen, maxlen)
        list = []
        for p in prefix_sorted
            for s in suffix_sorted
                push!(list, (prefix=p[1],prefix_count=p[2],suffix=s[1],suffix_count=s[2]))
            end
        end
        flank_sorted = sort(DataFrame(list),[:prefix_count,:suffix_count],rev=true)
        return map(x->Regex(".*(?<=$(x.prefix))(.{$minlen,$maxlen}?)(?=$(x.suffix))"), eachrow(flank_sorted))
    end

    """
        first_match(flank_regex, read)
    
    Find first match
    """
    function first_match(flank_regex, read)
        for r in flank_regex
            m = match(r, read)
            if m !== nothing
                return m
            end
        end
        return nothing
    end

    """
    Use regex to extract sequences between flanks
    """
    function search_regex(table, plist, slist; trim_ends=30, flank=100, minlen=15, maxlen=40)
        p = Progress(nrow(table), desc = "Searching sequences with regex")
        flank_regex = create_regex(plist, slist, minlen, maxlen)
        result = Folds.map(1:nrow(table)) do i
            seq = table[i,:genomic_sequence]
            seq = seq[trim_ends+1:end-trim_ends] # Trim potential primer
            m = first_match(flank_regex, seq)
            if m !== nothing
                hit = m[1]
                start = minimum(findfirst(m[1], seq))
                stop = maximum(findfirst(m[1], seq))
                prefix = seq[maximum([1,start-flank]):start-1]
                suffix = seq[stop+1:minimum([stop+flank,length(seq)])]
                return (row=i, hit=hit, prefix=prefix, suffix=suffix)
            else
                return (row=i, hit=nothing, prefix=nothing, suffix=nothing)
            end
            next!(p)
        end
        finish!(p)
        return result
    end


    """
        local_align(x, y)

    Align two sequences locally
    """
    function local_align(x, y)
        scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-5, gap_extend=-1);
        res = pairalign(SemiGlobalAlignment(), x, y, scoremodel)
    end
    
    """
        best_alignment(query, db)
    
    Find best alignment between query and database
    """
    function best_alignment(query, db)
        best_score = 0
        best_name = nothing
        best_aln = nothing
        best_ref = nothing
        best_distance = Inf
        best_regex_seq = nothing
        best_gaps = 0
        for (name, seq) in db
            aln = local_align(query, seq)
            s = score(aln)
            aln_str = last(split(string(aln.aln.a)))
            ref_str = last(split(string(aln.aln.b)))
            gaps = sum([x=='-' for x in aln_str])
            distance = sum([x!=y for (x,y) in zip(aln_str, ref_str)])
            if s > best_score
                best_score = s
                best_name = name
                best_ref = replace(ref_str,"-"=>"")
                best_aln = replace(aln_str,"-"=>"")
                best_regex_seq = query
                best_distance = distance
                best_gaps = gaps
            end
        end
        return (best_name = best_name, distance=best_distance, best_seq_gaps = best_gaps, best_score = best_score, regex_seq = best_regex_seq, best_aln = best_aln, best_ref = best_ref)
    end

    """
        annotate(agumented_table, db; nprefix=7, nsuffix=7)

    Annotate sequences with closest alleles
    """
    function annotate(agumented_table, db; nprefix=7, nsuffix=7)
        p = Progress(nrow(agumented_table), desc = "Alinging extracted sequences")
        result = Folds.map(eachrow(agumented_table)) do row
            best_name, distance, best_seq_gaps, best_score, regex_seq, best_aln, best_ref = best_alignment(row.hit, db)
            next!(p)
            return (well=row.well, case=row.case, best_name=best_name, distance=distance, gaps=best_seq_gaps, regex_seq=regex_seq, prefix=row.hit_prefix[end-nprefix+1:end], best_aln=best_aln, suffix=row.hit_suffix[1:nsuffix], best_ref=best_ref)
        end
        finish!(p)
        return DataFrame(result)
    end
end
