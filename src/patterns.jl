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
        kmers = [Set(split_kmers(group.seq[i], fragment_size)) for i in 1:length(group.seq)]
        common_kmers = intersect(kmers...)
        # Filter common_kmers if they occur within outgroup seq
        fragments = filter(x->!any([occursin(x, seq) for seq in outgroup.seq]), common_kmers)

        # Don't allow too long fragments as these might miss variation
        if isempty(fragments) && fragment_size < max_fragment_size
            @info "Trying fragment size $(fragment_size+1)"
            return find_unique_fragments(group::DataFrame, outgroup::DataFrame; fragment_size=fragment_size+1)
        end
        return unique(fragments)
    end

    function search_pattern(reads_df::DataFrame, group::DataFrame, outgroup::DataFrame, db_df::DataFrame, blacklist::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01, min_count=10, max_dist=10) where F <: AbstractFloat
        # Lookup dictionary for database sequences
        db_dict = Dict(map(x->(x[1],length(x[2])), eachrow(db_df)))
        # Find unique fragments
        fragments = find_unique_fragments(group, outgroup, fragment_size=fragment_size, max_fragment_size=max_fragment_size)
        if length(fragments) == 0
            @info "No patterns found for $(unique(group.gene))"
            return nothing
        end
        # Sample fragments
        search_fragments = unique(rand(fragments, max_fragments))

        # Filter df for rows which have a sequence that is a substring of any of the fragments
        gene_df = filter(x->any([occursin(fragment, x.genomic_sequence) for fragment in search_fragments]), reads_df)

        # Trim sequences
        prof_start, prof_stop = trim.trim_profiles(group.seq, weights)
        ok, trimmed_sequences, region = trim.find_segments(gene_df.genomic_sequence, prof_start, prof_stop, mlen)
        trimmed_sequences_df = DataFrame(seq = trimmed_sequences)
        clustered = sort(combine(groupby(trimmed_sequences_df, :seq), nrow => :count), :count, rev=true)

        # Annotate closest
        closest = map(x->find_closest(x, db_df), clustered.seq)
        clustered[!,:closest] = map(x->x[1], closest)
        clustered[!,:dist] = map(x->x[2], closest)
        clustered[!,:allele_name] = map(x->(x.dist == 0) ? x.closest : hamming.unique_name(x.closest,x.seq)*" Novel", eachrow(clustered))
        clustered[!,:gene] = map(x->first(split(x.allele_name,'*',limit=2)), eachrow(clustered))
        clustered[!,:length] = length.(clustered.seq)
        clustered[!,:length_db] = map(x->get(db_dict, x[1], "") , closest)
        clustered[!,:patterns] .= join(search_fragments, ",")
        # Compute and filter by frequency
        clustered[!,:ratio] = clustered.count ./ maximum(clustered.count)
        return clustered
    end

    function filter_and_report!(table::DataFrame, blacklist::DataFrame)
        # Define a local function that checks if any sequence in blacklist occurs in a given element's sequence
        function should_filter(x)
            if any(occursin.(x.seq, blacklist.seq))
                println("Discarding blacklisted: $(x.allele_name) with count $(x.count): $(x.seq)")
                return true   
            else
                return false
            end
        end
    
        # Use the filter! function with the local should_filter function
        filter!(x -> !should_filter(x), table)
    end
    
    # Write a function that runs search_pattern for each gene
    function search_patterns(reads_df::DataFrame, blacklist::DataFrame, db_df::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01, min_count=10, max_dist=10) where F <: AbstractFloat
        candidates = Vector{DataFrame}()
        for gene in unique(db_df.gene)
            @info gene
            group = filter(x->x.gene == gene, db_df)
            local outgroup
            if nrow(blacklist) > 0
                outgroup = reduce(vcat,[filter(x->x.gene != gene, db_df), filter(x->x.gene != gene, blacklist)])
            else
                outgroup = filter(x->x.gene != gene, db_df)
            end
            discovered = search_pattern(reads_df, group, outgroup, db_df, blacklist, fragment_size=fragment_size, max_fragment_size=max_fragment_size, max_fragments=max_fragments, weights=weights, mlen=mlen, min_freq=min_freq, min_count=min_count, max_dist=max_dist)
            if discovered !== nothing
                push!(candidates, discovered)
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
