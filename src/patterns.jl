module patterns
    include("trim.jl")
    include("data.jl")
    using .trim
    using .data
    using DataFrames
    using CSV
    using StringDistances
    using Folds
    using Logging
    using ProgressMeter

    export split_kmers, find_closest, find_unique_fragments, trim_by_distance, search_pattern, search_patterns

    # Write function splitting string into kmers
    function split_kmers(s::AbstractString, k::Int)
        kmers = Vector{String}()
        for i in 1:length(s)-k+1
            push!(kmers, s[i:i+k-1])
        end
        # Filter out kmers that appear more than once in a sequence
        unique_kmers = filter(y->count(x->x == y, kmers) == 1, kmers)
        return unique_kmers
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
        if length(kmers) == 0
            return []
        end
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

    """
        trim_by_distance(group::DataFrame, gene_df::DataFrame, search_fragments::Vector{String})

    Trim sequences by distance to 5' and 3' end and position of the search_fragments.
    The distance to 5' and 3' is determine from gene lengths in the group dataframe.
    """
    function trim_by_distance(group::DataFrame, gene_df::DataFrame, search_fragments::Vector{String})
        # Find positions of search_fragments in group dataframe
        positions = Dict()
        trimmed_sequences = []
        for kmer in search_fragments
            positions[kmer] = ()
        end
        for i in 1:nrow(group)
            for kmer in search_fragments
                pos = findfirst(kmer, group[i,:seq])
                if pos === nothing
                    continue
                end
                # Compute number of characters before and after the search_fragment
                prefix = minimum(pos)-1
                suffix = length(group[i,:seq]) - prefix - length(kmer)
                positions[kmer] = (prefix, suffix)
            end
        end
        # Find positions of search_fragments in gene_df dataframe
        for i in 1:nrow(gene_df)
            for kmer in search_fragments
                pos = findfirst(kmer, gene_df[i,:genomic_sequence])
                read_length = length(gene_df[i,:genomic_sequence])
                if pos === nothing
                    continue
                end
                prefix, suffix = positions[kmer]
                start = minimum(pos)
                if (start - prefix < 1) || (start + length(kmer) + suffix > read_length)
                    continue
                end
                gene = gene_df[i,:genomic_sequence][start-prefix:start+suffix+length(kmer)-1]
                push!(trimmed_sequences, gene)
            end
        end
        return trimmed_sequences
    end

    function search_pattern(reads_df::DataFrame, group::DataFrame, outgroup::DataFrame, db_df::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01, min_count=10, max_dist=10, noprofile=false) where F <: AbstractFloat
        # Lookup dictionary for database sequences
        db_dict = Dict(map(x->(x[1],length(x[2])), eachrow(db_df)))
        # Find unique fragments
        fragments = find_unique_fragments(group, outgroup, fragment_size=fragment_size, max_fragment_size=max_fragment_size)
        if length(fragments) == 0
            return nothing
        end
        # Sample fragments
        search_fragments = unique(rand(fragments, max_fragments))

        # Filter df for rows which have a sequence that is a substring of any of the fragments
        gene_df = filter(x->any([occursin(fragment, x.genomic_sequence) for fragment in search_fragments]), reads_df)
        attempt = 1
        while (nrow(gene_df) <= min_count)
            if attempt == 3
                @info "No patterns found for $(unique(group.gene))"
                return nothing
            end
            @info "Attempt $attempt at patterns resampling because no sufficient (>= $min_count) matches for $(unique(group.gene))"
            attempt += 1
            search_fragments = unique(rand(fragments, max_fragments))
            gene_df = filter(x->any([occursin(fragment, x.genomic_sequence) for fragment in search_fragments]), reads_df)
        end
        trimmed_sequences = []
        if !noprofile
            # Trim sequences
            prof_start, prof_stop = trim.trim_profiles(group.seq, weights)
            ok, trimmed_sequences, region = trim.find_segments(gene_df.genomic_sequence, prof_start, prof_stop, mlen)
            if sum(ok)/length(ok) < 0.9
                @warn "Not enough sequences trimmed for $(unique(group.gene)) possibly due incomplete database?"
            end
        else
            trimmed_sequences = trim_by_distance(group, gene_df, search_fragments)
        end
        # Group by identity
        trimmed_sequences_df = DataFrame(seq = trimmed_sequences)
        clustered = sort(combine(groupby(trimmed_sequences_df, :seq), nrow => :count), :count, rev=true)

        # Annotate closest
        closest = map(x->find_closest(x, db_df), clustered.seq)
        clustered[!,:closest] = map(x->x[1], closest)
        clustered[!,:dist] = map(x->x[2], closest)
        clustered[!,:allele_name] = map(x->(x.dist == 0) ? x.closest : unique_name(x.closest,x.seq)*" Novel", eachrow(clustered))
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
            # Check if seq contains any characters that are ambigous
            if any([c ∉ ['A','T','G','C'] for c in x.seq])
                println("Discarding sequence with ambigous nucleotides: $(x.allele_name) with count $(x.count): $(x.seq)")
                return true
            end
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
    function search_patterns(reads_df::DataFrame, blacklist::DataFrame, db_df::DataFrame; fragment_size::Int=20, max_fragment_size::Int=60, max_fragments::Int=5, weights::Int=20, mlen::Int=100, min_freq::F=0.01, min_count=10, max_dist=10, noprofile=false) where F <: AbstractFloat
        # Assert that group, outgroup and db_df table have required colums
        @assert all([x in names(blacklist) for x in ["id", "gene", "seq"]])
        @assert all([x in names(db_df) for x in ["id", "gene", "seq"]])

        # Initialize an array to hold loggers for each thread
        thread_loggers = [Logging.SimpleLogger(IOBuffer()) for _ in 1:Threads.nthreads()]
        p = ProgressMeter.Progress(length(unique(db_df.gene)), dt=0.5, showspeed=true)
        candidates = Folds.map(unique(db_df.gene)) do gene
            logger = thread_loggers[Threads.threadid()]
            with_logger(logger) do
                @info gene
                group = filter(x->x.gene == gene, db_df)
                local outgroup
                if nrow(blacklist) > 0
                    outgroup = reduce(vcat,[filter(x->x.gene != gene, db_df), filter(x->x.gene != gene, blacklist)])
                else
                    outgroup = filter(x->x.gene != gene, db_df)
                end
                discovered = search_pattern(reads_df, group, outgroup, db_df, fragment_size=fragment_size, max_fragment_size=max_fragment_size, max_fragments=max_fragments, weights=weights, mlen=mlen, min_freq=min_freq, min_count=min_count, max_dist=max_dist, noprofile=noprofile)
                ProgressMeter.next!(p)
                if discovered !== nothing
                    return discovered
                else
                    @warn "No patterns found for $gene"
                    return DataFrame([[],[],[],[], [],[],[], [],[],[]], [:count, :closest, :dist, :allele_name, :gene, :length, :length_db, :patterns, :ratio, :seq])
                end
            end
        end
        ProgressMeter.finish!(p)
        # Combine and display the log outputs at the end
        for logger in thread_loggers
            seekstart(logger.stream)
            print(read(logger.stream, String))
        end
        final = reduce(vcat, candidates)
        filter_and_report!(final, blacklist)
        other_columns = setdiff(names(final), ["seq"])
        final = final[:, vcat(other_columns, "seq")]
        return final
    end
end
