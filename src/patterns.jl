module patterns
    include("trim.jl")
    include("data.jl")
    include("hamming.jl")
    using .trim
    using .data
    using .hamming
    using DataFrames
    using CSV
    using StringDistances

    # Write function splitting string into kmers
    function split_kmers(s::AbstractString, k::Int)
        kmers = Vector{String}()
        for i in 1:length(s)-k+1
            push!(kmers, s[i:i+k-1])
        end
        return kmers
    end

    function find_closest(seq::AbstractString, db_df::DataFrame)
        min_dist = 100
        closest = ""
        edit_dist = Levenshtein()
        for i in 1:nrow(db_df)
            dist = edit_dist(seq, db_df[i,:seq])
            if dist < min_dist
                min_dist = dist
                closest = db_df[i,:id]
            end
        end
        return (closest, min_dist)
    end

    function find_unique_fragments(group::DataFrame, outgroup::DataFrame; fragment_size::Int=15, max_fragment_size::Int=60)
        # Search for fragments
        @info "Group length $(length(group.seq))"
        kmers = [Set(split_kmers(group.seq[i], fragment_size)) for i in 1:length(group.seq)]
        common_kmers = intersect(kmers...)
        # Filter common_kmers if they occur within outgroup seq
        fragments = filter(x->!any([occursin(x, seq) for seq in outgroup.seq]), common_kmers)

        # Don't allow too long fragments as these might miss variation
        if isempty(fragments) && fragment_size <= max_fragment_size
            return find_unique_fragments(group::DataFrame, outgroup::DataFrame; fragment_size=fragment_size+1)
        end
        return unique(fragments)
    end

    function search_pattern(df::DataFrame, group::DataFrame, outgroup::DataFrame, db_df::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01) where F <: AbstractFloat
        # Find unique fragments
        fragments = find_unique_fragments(group, outgroup, fragment_size=fragment_size, max_fragment_size=max_fragment_size)
        if length(fragments) == 0
            return nothing
        end
        # Sample fragments
        search_fragments = unique(rand(fragments, max_fragments))

        # Filter df for rows which have a sequence that is a substring of any of the fragments
        gene_df = filter(x->any([occursin(fragment, x.genomic_sequence) for fragment in search_fragments]), df)
    
        # Trim sequences
        prof_start, prof_stop = trim.trim_profiles(group.seq, weights)
        ok, trimmed_sequences, region = trim.find_segments(gene_df.genomic_sequence, prof_start, prof_stop, mlen)
        trimmed_sequences_df = DataFrame(seq = trimmed_sequences)
        clustered = sort(combine(groupby(trimmed_sequences_df, :seq), nrow => :count), :count, rev=true)
    
        # Compute and filter by frequency
        clustered[!,:frequency] = clustered.count ./ sum(clustered.count)
        clustered = filter(x->x.frequency >= min_freq, clustered)
        
        # Annotate closest
        closest = map(x->find_closest(x, db_df), clustered.seq)
        clustered[!,:closest] = map(x->x[1], closest)
        clustered[!,:dist] = map(x->x[2], closest)
        clustered[!,:allele_name] = map(x->(x.dist == 0) ? x.closest : hamming.unique_name(x.closest,x.seq)*" Novel", eachrow(clustered))
        clustered[!,:patterns] .= join(search_fragments, ",")
        return clustered
    end

    # Write a function that runs search_pattern for each gene
    function search_patterns(df::DataFrame, blacklist::DataFrame, db_df::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01) where F <: AbstractFloat
        candidates = Vector{DataFrame}()
        for gene in unique(db_df.gene)
            @info gene
            group = filter(x->x.gene == gene, db_df)
            local outgroup
            if nrow(blacklist) > 0
                @info "Using blacklist"
                outgroup = reduce(vcat,[filter(x->x.gene != gene, group), filter(x->x.gene != gene, blacklist)])
            else
                @info "Not using blacklist"
                outgroup = filter(x->x.gene != gene, group)
                @info "Outgroup $(nrow(outgroup))"
            end
            found_patterns = search_pattern(df, group, outgroup, db_df, fragment_size=fragment_size, max_fragment_size=max_fragment_size, max_fragments=max_fragments, weights=weights, mlen=mlen, min_freq=min_freq)
            if found_patterns !== nothing
                push!(candidates, found_patterns)
            else
                @warn "No patterns found for $gene"
            end
        end
        final = reduce(vcat, candidates)
        other_columns = setdiff(names(final), ["seq"])
        final = final[:, vcat(other_columns, "seq")]
        return final
    end
end
