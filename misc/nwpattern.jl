module nwpattern
    using DataFrames
    using CSV
    using StringDistances
    using Folds
    using Logging
    using ProgressMeter
    using FASTX
    using BioSequences
    using Kmers
    using DataStructures
    using Combinatorics
    using Statistics
    using Random
    using Folds
    include("trim.jl")
    include("data.jl")
    export discover_alleles, load_and_filter_csv, load_fasta_biosequences, create_hashtable, find_indices_for_priority_kmers, find_gene_unique_kmers, trim_sequences, find_closest

    """
        load_fasta_biosequences(path::String)::Vector{Tuple{String,BioSequence}}

    Load the FASTA file as DNA sequences
    """
    function load_fasta_biosequences(path::String)::Vector{Tuple{String,BioSequence}}
        records = Vector{Tuple{String,BioSequence}}()
        open(FASTA.Reader, path) do reader
            for record in reader
                push!(records, (FASTA.identifier(record), FASTA.sequence(LongDNA{2}, record)))
            end
        end
        return records
    end

    """
        create_hashtable(sequences, k)

    Create a hashtable of kmers from the sequences
    """
    function create_hashtable(sequences, k)
        hashtable = Dict{DNAKmer{k}, BitSet}()
        len = length(sequences)
        n = Progress(len, desc="Creating hashtable for $len sequences")
        @inbounds for (i, seq) in enumerate(sequences)
            @inbounds for (j, kmer) in EveryKmer{DNAKmer{k}}(seq)
                set = get!(hashtable, kmer, BitSet())
                push!(set, i)
            end
            next!(n)
        end
        finish!(n)
        return hashtable
    end

    """
        load_and_filter_csv(path::String; delim::Char='\t')

    Load and filter the CSV file
    """
    function load_and_filter_csv(path::String; delim::Char='\t')
        data = CSV.File(path, delim=delim) |> DataFrame
        filter!(x -> !occursin('N', x.genomic_sequence) && length(x.genomic_sequence) != 0, data)
        return data
    end

    """
        find_indices_for_priority_kmers(priority, genomic_hashmap; n=3, max_combinations=100)

    Find the combination of n kmers such that the combination is not present in the outgroup priority list
    """
    function find_indices_for_priority_kmers(priority, genomic_hashmap; n=3, max_combinations=100)
        indices = Set{Int}()
        combi = shuffle(collect(combinations(first(priority,100), n)))

        for comb in first(combi, max_combinations)
            for i in intersect([get(genomic_hashmap, c, "") for c in comb]...)
                push!(indices, i)
            end
        end
    end

    function find_gene_unique_kmers(kmers, genomic_hashmap, n; max_combinations=10000, max_attempts = 10000)
        indices = Set{Int}()
        used_combinations = Set()
        attempt = 0
        while length(used_combinations) < max_combinations
            comb = Set(first(shuffle(kmers), n))
            attempt += 1
            if (length(used_combinations) >= max_combinations) | (attempt >= max_attempts)
                break
            end
            # Retry if combination already used
            if comb in used_combinations
                continue
            end
            if length(comb) != n
                continue
            end
            # Check for combination
            for i in intersect([get(genomic_hashmap, c, "") for c in comb]...)
                push!(indices, i)
                push!(used_combinations, comb)
            end
        end
        
        return indices, used_combinations
    end

    """
        trim_sequences(gene_data, target; flank_length=21, min_length=200)

    Trim the sequences based on the target sequences
    """
    function trim_sequences(gene_data, target; flank_length=21, min_length=200)
        prof_start, prof_stop = trim.trim_profiles(target, flank_length)
        ok, trimmed_sequences, region = trim.find_segments(gene_data.genomic_sequence, prof_start, prof_stop, min_length)
        gene_data = gene_data[ok,:]
        gene_data[:, :sequence] = trimmed_sequences
        return gene_data
    end

    """
        find_closest(seq::AbstractString, db_df::DataFrame)

    Find the closest sequence to the given sequence
    """
    function find_closest(seq::AbstractString, db::AbstractVector)
        min_dist = 100
        closest = ""
        edit_dist = Levenshtein()
        for i in eachindex(db)
            dist = edit_dist(seq, string(last(db[i])))
            if dist < min_dist
                min_dist = dist
                closest = first(db[i])
            end
        end
        return (closest, min_dist)
    end

    function discover_alleles(tsv_path, fasta_path, candidates_tsv, filtered_candidates_tsv; k=7, n=3, max_distance=10, flank_length=21, min_length=200, min_count=5, min_frequency=0.1, max_combinations=10000, max_attempts=10000)
        table = load_and_filter_csv(tsv_path)
        db = load_fasta_biosequences(fasta_path)
        db_genes = map(x->first(split(first(x),'*')), db)
        db_lookup = Dict([(n, string(s)) for (n,s) in db])
        genomic_hashmap = create_hashtable(map(LongDNA{2}, table.genomic_sequence), k)
        candidates = []
        for gene in unique(db_genes)
            @info "------------- Gene $gene -------------"
            start_time = time()
            db_gene = db[db_genes .== gene]
            group_hash = create_hashtable(map(x->last(x), db_gene), k)
            outgroup_hash = create_hashtable(map(x->last(x), db[.!(db_genes .== gene)]), k)
            kmers = collect(setdiff(keys(group_hash), keys(outgroup_hash)))
            indices, used_comb = find_gene_unique_kmers(kmers, genomic_hashmap, n, max_combinations=max_combinations, max_attempts=max_attempts)
            gene_data = trim_sequences(table[collect(indices),:], string.(last.(db_gene)), flank_length=flank_length, min_length=min_length)
            # Assign gene and count reads per well and case
            gene_data_df = sort(combine(groupby(gene_data, [:well, :case, :sequence]), nrow => :count), :count, rev=true)
            end_time = time()
            elapsed_time = end_time - start_time
            if length(kmers) == 0
                @warn "Skipping gene: $gene, no gene specific kmers found!"
                continue
            elseif length(indices) == 0
                @warn "Skipping gene: $gene, no gene specific reads matching kmer combination found!"
                continue
            end
            # NOTE: indices are integers not boolean values
            @info "Gene: $gene, Unique kmers: $(length(kmers)), kmer combinations: $(length(used_comb)), matched sequences $(length(indices)), max_case_count $(maximum(gene_data_df.count)), took $(round(elapsed_time, digits=3)) seconds."
            gene_data_df[:,:gene] .= gene
            prog = Progress(nrow(gene_data_df), desc="Finding closest sequences for gene: $gene")
            closest_and_distances = Folds.map(eachrow(gene_data_df)) do row
                ProgressMeter.next!(prog)
                find_closest(row.sequence, db_gene)
            end
            ProgressMeter.finish!(prog)
            gene_data_df[!, :closest] = first.(closest_and_distances)
            gene_data_df[!, :distance] = last.(closest_and_distances)
    
            # Add column containing sequence length
            transform!(gene_data_df, :sequence => ByRow(length) => :length)
            # Add column containing closest sequence length
            transform!(gene_data_df, :closest => ByRow(x -> length(get(db_lookup,x,""))) => :closest_length)
            # Add number of sequences matching all kmer combinations
            transform!(gene_data_df, :sequence => ByRow(x -> length(indices)) => :gene_matches)
            # Add allelic ratio frequency
            transform!(groupby(gene_data_df, [:well, :case, :gene]), :count => (x -> x ./ maximum(x)) => :frequency)
            
            # Collect the candidates
            push!(candidates, gene_data_df)
        end

        # Combine all genes
        candidates_df = reduce(vcat, candidates)

        # Compute new name for all rows where distance != 0 based on column closest
        transform!(candidates_df, [:distance, :closest, :sequence] => 
        ((distance, closest, sequence) -> begin
            map((d, c, s) -> d == 0 ? c : data.unique_name(c, s) * " Novel", distance, closest, sequence)
        end) => :allele_name)
        
        candidates_df_all = sort(candidates_df, [:well, :case, :gene,:count], rev=[false, false, false, true])
        CSV.write(candidates_tsv, candidates_df_all, delim='\t', compress=true)
        @info "Saved candidates to $candidates_tsv file."
    
        # Filter the candidates based on min_count and min_frequency
        candidates_df_filtlered = filter(x -> (x.count >= min_count) && (x.frequency >= min_frequency) && (x.distance <= max_distance), candidates_df)
        CSV.write(filtered_candidates_tsv, candidates_df_filtlered, delim='\t', compress=true)
        @info "Saved filtered candidates to $filtered_candidates_tsv file."
    end
end
