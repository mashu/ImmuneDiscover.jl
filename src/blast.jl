module blast
    using CSV
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    using BioSequences
    using Folds

    using ..data: load_fasta as data_load_fasta, unique_name

    export blast_discover, save_to_fasta, accumulate_affixes, save_extended, handle_blast

    const columns = ["qseqid", "sseqid", "pident", "nident", "length", "mismatch", "gapopen", "qcovs", "qcovhsp", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore", "sstrand", "qseq"]

    # --- AlignmentStats tracks prefix/suffix trimming outcomes ---

    mutable struct AlignmentStats
        total_attempts::Int
        prefix_failures::Int
        suffix_failures::Int
        successful_trims::DefaultDict{String, Vector{String}}
        failed_genes::DefaultDict{String, Vector{String}}
    end

    AlignmentStats() = AlignmentStats(0, 0, 0,
        DefaultDict{String, Vector{String}}(Vector{String}),
        DefaultDict{String, Vector{String}}(Vector{String}))

    """
        read_fasta(filepath) -> Vector{Tuple{String, String}}

    Delegates to data.load_fasta for shared FASTA-reading logic.
    """
    function read_fasta(filepath::String)
        return data_load_fasta(filepath, validate=false)
    end

    load_csv(path::String; delim::Char='\t') = CSV.File(path, delim=delim) |> DataFrame
    select_columns(df::DataFrame, cols::Vector{Symbol}) = unique(select(df, cols))

    function save_to_fasta(records, output_file::String)
        open(output_file, "w") do io
            for (well, case, name, sequence) in records
                write(io, ">$name $well $case\n$sequence\n")
            end
        end
    end

    function save_to_fasta(records::Vector{Tuple{String, String}}, output_file::String)
        open(output_file, "w") do io
            for (name, sequence) in records
                write(io, ">$name \n$sequence\n")
            end
        end
    end

    function time_command(cmd::Cmd)
        start_time = time()
        run(cmd)
        return time() - start_time
    end

    function blastn(query_file::String, database::String, output_file::String; args::String="")
        outfmt = "6 " * join(columns, " ")
        cmd = `blastn -query $query_file -subject $database -out $output_file -outfmt $outfmt`
        if !isempty(args)
            cmd = `$cmd $(split(args))`
        end
        elapsed = time_command(cmd)
        @info "BLASTn completed in $(round(elapsed, digits=2)) seconds"
    end

    function replace_extension(path, ext)
        dir = dirname(path)
        base = first(split(basename(path), '.'))
        return joinpath(dir, base * "." * ext)
    end

    function longest_common_suffix_str(a::String, b::String)
        i = 0
        while i < min(length(a), length(b)) && a[end-i] == b[end-i]
            i += 1
        end
        return a[end-i+1:end]
    end

    function longest_common_prefix_str(a::String, b::String)
        i = 0
        while i < min(length(a), length(b)) && a[i+1] == b[i+1]
            i += 1
        end
        return a[1:i]
    end

    function accumulate_affixes(db, demux_df; forward_extension=20, reverse_extension=20)
        singleton = Vector{Tuple{String, String, String, String}}()
        for (name, reference_seq) in db
            matches = filter(x -> occursin(reference_seq, x.genomic_sequence), eachrow(demux_df))
            if length(matches) == 0
                push!(singleton, (name, reference_seq, "", ""))
                continue
            end
            common_prefix = ""
            common_suffix = ""
            for row in matches
                gs = row.genomic_sequence
                m = findfirst(reference_seq, gs)
                m === nothing && continue
                start_pos = minimum(m)
                end_pos = maximum(m)
                prefix = gs[max(1, start_pos - forward_extension):start_pos - 1]
                suffix = gs[end_pos + 1:min(end_pos + reverse_extension, length(gs))]
                if common_prefix == ""
                    common_prefix = prefix
                    common_suffix = suffix
                else
                    common_prefix = longest_common_suffix_str(common_prefix, prefix)
                    common_suffix = longest_common_prefix_str(common_suffix, suffix)
                end
            end
            extended_sequence = common_prefix * reference_seq * common_suffix
            push!(singleton, (name, extended_sequence, common_prefix, common_suffix))
        end
        return singleton
    end

    function save_extended(extended_Ds, fasta_path)
        base_affixes = Vector{Tuple{String, String, String}}()
        open(fasta_path, "w") do io
            for (name, sequence, prefix, suffix) in extended_Ds
                write(io, ">$name\n$sequence\n")
                push!(base_affixes, (name, prefix, suffix))
            end
        end
        return base_affixes
    end

    nogaps(s) = replace(s, '-' => "")

    function check_affix_quality_warning(affix_length::Int, quality_threshold::Float64)
        if affix_length > 0 && affix_length <= 20 && quality_threshold > 0.5
            @warn "Quality threshold $(round(quality_threshold * 100, digits=1))% might be too strict for short affixes ($affix_length nt). Consider lowering --minquality."
        end
    end

    function calculate_alignment_quality(pairs, affix_is_first::Bool=true)
        if affix_is_first
            affix_start = findfirst(p -> first(p) != DNA_Gap, pairs)
            affix_end = findlast(p -> first(p) != DNA_Gap, pairs)
        else
            affix_start = findfirst(p -> last(p) != DNA_Gap, pairs)
            affix_end = findlast(p -> last(p) != DNA_Gap, pairs)
        end
        (affix_start === nothing || affix_end === nothing) && return 0.0
        aligned_region = pairs[affix_start:affix_end]
        matches = sum(first.(aligned_region) .== last.(aligned_region))
        return matches / length(aligned_region)
    end

    """
        safe_pairalign(seq1, seq2, scoremodel)

    Align two sequences, returning nothing if either is empty.
    """
    function safe_pairalign(seq1::LongDNA{4}, seq2::LongDNA{4}, scoremodel::AffineGapScoreModel)
        if length(seq1) == 0 || length(seq2) == 0
            @warn "Alignment skipped: one or both sequences are empty"
            return nothing
        end
        return pairalign(SemiGlobalAlignment(), seq1, seq2, scoremodel)
    end

    function remove_gaps(query::LongDNA{4})::LongDNA{4}
        return LongDNA{4}(filter(nt -> nt != DNA_Gap, query))
    end

    function trim_sequence(query::LongDNA{4}, prefix::LongDNA{4}, suffix::LongDNA{4}, stats,
        scoremodel::AffineGapScoreModel=AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1);
        min_quality=0.75, sseqid="")

        if length(query) == 0
            push!(stats.failed_genes[sseqid], "Empty query sequence")
            return nothing
        end
        partial_query = query
        stats.total_attempts += 1

        if length(prefix) > 0
            prefix_aln_result = safe_pairalign(prefix, query, scoremodel)
            if prefix_aln_result === nothing
                push!(stats.failed_genes[sseqid], "Prefix alignment failed")
                stats.prefix_failures += 1
                return nothing
            end
            prefix_aln = alignment(prefix_aln_result)
            prefix_pairs = collect(prefix_aln)
            prefix_positions = findall(p -> first(p) != DNA_Gap, prefix_pairs)
            if isempty(prefix_positions)
                push!(stats.failed_genes[sseqid], "No prefix content in alignment")
                stats.prefix_failures += 1
                return nothing
            end
            matches = sum(first(prefix_pairs[i]) == last(prefix_pairs[i]) for i in prefix_positions)
            prefix_quality = matches / length(prefix_positions)
            if prefix_quality < min_quality
                push!(stats.failed_genes[sseqid], "Poor prefix alignment quality ($(round(prefix_quality * 100, digits=1))% match)")
                stats.prefix_failures += 1
                return nothing
            end
            aligned_prefix = first.(prefix_pairs)
            aligned_query = last.(prefix_pairs)
            prefix_end = findlast(x -> x != DNA_Gap, aligned_prefix)
            if prefix_end === nothing
                push!(stats.failed_genes[sseqid], "No prefix match found")
                stats.prefix_failures += 1
                return nothing
            end
            query_core_start = prefix_end + 1
            if query_core_start > length(aligned_query)
                push!(stats.failed_genes[sseqid], "Query too short after prefix")
                stats.prefix_failures += 1
                return nothing
            end
            core_part = aligned_query[query_core_start:end]
            partial_query = remove_gaps(LongDNA{4}(join(core_part)))
            if length(partial_query) == 0
                push!(stats.failed_genes[sseqid], "Empty core after prefix trimming")
                stats.prefix_failures += 1
                return nothing
            end
        end

        if length(suffix) > 0
            suffix_aln_result = safe_pairalign(suffix, partial_query, scoremodel)
            if suffix_aln_result === nothing
                push!(stats.failed_genes[sseqid], "Suffix alignment failed")
                stats.suffix_failures += 1
                return nothing
            end
            suffix_aln = alignment(suffix_aln_result)
            suffix_pairs = collect(suffix_aln)
            suffix_positions = findall(p -> first(p) != DNA_Gap, suffix_pairs)
            if isempty(suffix_positions)
                push!(stats.failed_genes[sseqid], "No suffix content in alignment")
                stats.suffix_failures += 1
                return nothing
            end
            matches = sum(first(suffix_pairs[i]) == last(suffix_pairs[i]) for i in suffix_positions)
            suffix_quality = matches / length(suffix_positions)
            if suffix_quality < min_quality
                push!(stats.failed_genes[sseqid], "Poor suffix alignment quality ($(round(suffix_quality * 100, digits=1))% match)")
                stats.suffix_failures += 1
                return nothing
            end
            aligned_suffix = first.(suffix_pairs)
            aligned_query = last.(suffix_pairs)
            suffix_start = findfirst(x -> x != DNA_Gap, aligned_suffix)
            if suffix_start === nothing
                push!(stats.failed_genes[sseqid], "No suffix match found")
                stats.suffix_failures += 1
                return nothing
            end
            if suffix_start <= 1
                push!(stats.failed_genes[sseqid], "Suffix starts too early")
                stats.suffix_failures += 1
                return nothing
            end
            core_part = aligned_query[1:suffix_start-1]
            partial_query = remove_gaps(LongDNA{4}(join(core_part)))
            if length(partial_query) == 0
                push!(stats.failed_genes[sseqid], "Empty core after suffix trimming")
                stats.suffix_failures += 1
                return nothing
            end
        end

        core_length = length(partial_query)
        if core_length < 10
            push!(stats.failed_genes[sseqid], "Core sequence too short ($core_length < 10 nt)")
            stats.prefix_failures += 1
            return nothing
        end

        push!(stats.successful_trims[sseqid], "Trimmed successfully")
        return partial_query
    end

    function compute_edit_distance(query::String, reference::String)
        aln = pairalign(LevenshteinDistance(), reference, query)
        return score(aln)
    end

    function trim_and_align_sequence(query::String, prefix::String, suffix::String, reference::String, stats; min_quality=0.75, sseqid="")
        query_dna = BioSequences.LongDNA{4}(query)
        prefix_dna = BioSequences.LongDNA{4}(prefix)
        suffix_dna = BioSequences.LongDNA{4}(suffix)
        trimmed_dna = trim_sequence(query_dna, prefix_dna, suffix_dna, stats, min_quality=min_quality, sseqid=sseqid)
        if trimmed_dna === nothing
            return "", -1
        end
        trimmed = replace(String(trimmed_dna), '-' => "")
        distance = compute_edit_distance(trimmed, reference)
        return trimmed, distance
    end

    """
        verify_blastn_version(min_version, max_version)

    Check that blastn is available and within the supported version range.
    """
    function verify_blastn_version(min_version::VersionNumber, max_version::VersionNumber)
        blastn_path = Sys.which("blastn")
        if blastn_path === nothing
            return (false, "blastn command not found in PATH")
        end
        cmd_output = read(`blastn -version`, String)
        version_match = match(r"blastn:\s+(\d+\.\d+\.\d+)", cmd_output)
        if isnothing(version_match)
            return (false, "Could not parse blastn version")
        end
        current_version = VersionNumber(version_match[1])
        if min_version ≤ current_version ≤ max_version
            return (true, "blastn version $current_version is within acceptable range")
        else
            return (false, "blastn version $current_version is outside acceptable range ($min_version - $max_version)")
        end
    end

    function edge(qseq, read)
        m = findfirst(nogaps(qseq), read)
        isnothing(m) && return missing, missing
        five_prime, three_prime = extrema(m)
        return five_prime, length(read) - three_prime
    end

    """
        blast_discover(tsv_path, combined_db_fasta; kwargs...)

    Perform assignments and discovery of alleles based on BLAST results.
    """
    function blast_discover(tsv_path, combined_db_fasta; max_dist=10, min_edge=10, min_scov=0.1, args="", verbose=false, overwrite=false)
        min_ver = v"2.15.0"
        max_ver = v"2.17.0"
        is_valid, message = verify_blastn_version(min_ver, max_ver)
        @info message
        if !is_valid
            error("Please install BLAST version $min_ver - $max_ver")
        end

        df = load_csv(tsv_path)
        @info "Read $(nrow(df)) rows from $tsv_path before BLAST assignment"
        query_fasta = replace_extension(tsv_path, "fasta")
        blast_tsv = replace_extension(tsv_path, "blast")
        unique_df = unique(df, :name)
        @info "Unique reads: $(nrow(unique_df))"
        if nrow(df) != nrow(unique_df)
            @error "Duplicated names found in input. Using unique ones but this is likely user error!"
        end
        df = unique_df

        query_sequences = collect.(eachrow(select_columns(df, [:well, :case, :name, :genomic_sequence])))
        save_to_fasta(query_sequences, query_fasta)

        blast_file = blast_tsv * ".gz"
        file_exists = isfile(blast_file)
        @info "BLAST cache check: file='$blast_file', exists=$file_exists, overwrite=$overwrite"
        if file_exists && !overwrite
            @info "BLASTn results already exist $blast_file. Skipping BLASTn."
        else
            @info "Running BLASTn. This may take a while."
            blastn(query_fasta, combined_db_fasta, blast_tsv, args=args)
            run(`gzip -f $blast_tsv`)
        end

        rm(query_fasta)

        blast_df = CSV.File(blast_tsv * ".gz", delim='\t', header=columns) |> DataFrame
        @info "BLASTn results read from $(blast_tsv).gz: $(nrow(blast_df)) rows"
        if nrow(blast_df) == 0
            error("No BLASTn results found (wrong BLAST parameters?). Cannot proceed.")
        end

        blast_df = combine(groupby(blast_df, :qseqid), x -> first(sort(x, [:pident, :qcovhsp, :qcovs, :bitscore], rev=true)))
        leftjoin!(blast_df, df, on=:qseqid => :name)
        transform!(blast_df, [:qseq, :genomic_sequence] => ByRow(edge) => [:five_prime_edge, :three_prime_edge])

        @info "BLASTn results after best hits: $(nrow(blast_df)) rows"
        filter!(x -> x.five_prime_edge > min_edge && x.three_prime_edge > min_edge, blast_df)
        @info "After filtering edge > $min_edge: $(nrow(blast_df)) rows"

        transform!(blast_df, [:length, :slen] => ByRow((len, slen) -> len / slen) => :scov)
        filter!(x -> x.scov > min_scov, blast_df)
        @info "After filtering scov > $min_scov: $(nrow(blast_df)) rows"

        transform!(blast_df, :qseq => ByRow(x -> replace(x, "-" => "")) => :qseq)
        filter!(x -> !startswith(x.sseqid, "P"), blast_df)
        @info "After filtering pseudo genes: $(nrow(blast_df)) rows"
        verbose && CSV.write(blast_tsv * "-pseudo.tsv", blast_df)

        read_name = Dict([(r.name, (r.well, r.case)) for r in eachrow(df)])
        blast_df[:, :well] = [read_name[x.qseqid][1] for x in eachrow(blast_df)]
        blast_df[:, :case] = [read_name[x.qseqid][2] for x in eachrow(blast_df)]

        clusters = combine(groupby(blast_df, [:well, :case, :sseqid, :qseq, :mismatch]), :qseqid => length => :full_count)
        @info "Clusters after grouping: $(nrow(clusters)) rows"
        verbose && CSV.write(blast_tsv * "-clusters.tsv", clusters)

        filter!(x -> x.mismatch <= max_dist, clusters)
        @info "After filtering mismatches <= $max_dist: $(nrow(clusters)) rows"
        verbose && CSV.write(blast_tsv * "-clusters-mismatch.tsv", clusters)

        transform!(clusters, :sseqid => ByRow(x -> split(x, "*")[1]) => :gene)
        transform!(groupby(clusters, [:well, :case, :gene]), :full_count => (x -> x ./ maximum(x)) => :full_frequency)

        return sort(clusters, [:well, :case, :sseqid], rev=false)
    end

    """
        handle_blast(parsed_args, immunediscover_module, always_gz)

    Handle the blast search pipeline. Extracted from real_main for modularity.
    """
    function handle_blast(parsed_args, immunediscover_module, always_gz)
        using_cli = immunediscover_module.cli
        parsed_args = using_cli.apply_blast_presets!(parsed_args)
        gene = parsed_args["search"]["blast"]["gene"]
        if haskey(using_cli.BLAST_PRESETS, gene)
            @info "Initial $gene preset parameters:"
            for (param, value) in using_cli.BLAST_PRESETS[gene]
                @info "  --$param = $value"
            end
        end
        @info "Running with the following parameters:"
        for (param, value) in parsed_args["search"]["blast"]
            @info "  --$param = $value"
        end
        @info "Discovery with BLAST assignments"

        fasta_path = parsed_args["search"]["blast"]["fasta"]
        work_dir = joinpath(dirname(fasta_path), ".immunediscover")
        isdir(work_dir) || mkpath(work_dir)
        file_stem = split(basename(fasta_path), '.')[1]
        affixes_path = joinpath(work_dir, file_stem * ".affixes")
        DB = immunediscover_module.load_fasta(fasta_path, validate=false)

        verbose = parsed_args["search"]["blast"]["verbose"]
        overwrite = parsed_args["search"]["blast"]["overwrite"]
        minquality = parsed_args["search"]["blast"]["minquality"]
        forward_extension = parsed_args["search"]["blast"]["forward"]
        reverse_extension = parsed_args["search"]["blast"]["reverse"]
        keep_failed = parsed_args["search"]["blast"]["keep-failed"]
        min_corecov = parsed_args["search"]["blast"]["min-corecov"]

        if (forward_extension < 7) && (forward_extension > 0)
            @warn "Forward extension $forward_extension is short and may lead to false positives"
        end
        if (reverse_extension < 7) && (reverse_extension > 0)
            @warn "Reverse extension $reverse_extension is short and may lead to false positives"
        end

        # Process pseudo genes if provided
        db_p = Vector{Tuple{String, String}}()
        pseudo = parsed_args["search"]["blast"]["pseudo"]
        if !isempty(pseudo)
            for (name, seq) in read_fasta(pseudo)
                push!(db_p, ("P" * name, seq))
            end
        end
        for (name, seq) in read_fasta(fasta_path)
            push!(db_p, (name, seq))
        end
        combined_fasta_path = joinpath(work_dir, file_stem * "-combined.fasta")
        save_to_fasta(db_p, combined_fasta_path)

        # Handle sequence extension
        if forward_extension == 0 && reverse_extension == 0
            ext_fasta_path = combined_fasta_path
            @info "No sequence extension requested, using original sequences"
            empty_affixes = [(name="", prefix="", suffix="")]
            CSV.write(affixes_path, DataFrame(empty_affixes), delim='\t')
        else
            ext_fasta_path = joinpath(work_dir, file_stem * "-combined-extended.fasta")
            if !isfile(ext_fasta_path) || overwrite
                @info "Extending gene sequences by $forward_extension forward and $reverse_extension reverse nucleotides"
                demux = load_csv(parsed_args["search"]["blast"]["input"])
                extended = accumulate_affixes(DB, demux,
                    forward_extension=forward_extension,
                    reverse_extension=reverse_extension)
                affixes = save_extended(extended, ext_fasta_path)
                CSV.write(affixes_path, DataFrame(affixes, [:name, :prefix, :suffix]), delim='\t')
                @info "Saved affixes in $affixes_path"
            else
                @info "Using existing extended sequences from $ext_fasta_path"
            end
        end

        # Run BLAST discovery
        blast_clusters = blast_discover(
            parsed_args["search"]["blast"]["input"],
            ext_fasta_path,
            max_dist=parsed_args["search"]["blast"]["maxdist"],
            min_edge=parsed_args["search"]["blast"]["edge"],
            min_scov=parsed_args["search"]["blast"]["subjectcov"],
            args=parsed_args["search"]["blast"]["args"],
            verbose=verbose,
            overwrite=overwrite
        )

        # Trim extensions if applicable
        if forward_extension != 0 || reverse_extension != 0
            if isfile(affixes_path)
                @info "Loading affixes from $affixes_path"
                affixes = Tuple.(collect.(eachrow(CSV.File(affixes_path, delim='\t') |> DataFrame)))
                affix_dict = Dict{String, Tuple{String, String}}()
                for (name, prefix, suffix) in affixes
                    affix_dict[name] = (prefix, suffix)
                end

                stats = AlignmentStats()
                check_affix_quality_warning(forward_extension, minquality)
                check_affix_quality_warning(reverse_extension, minquality)

                results = Folds.map(eachrow(blast_clusters)) do row
                    sseqid = String(row.sseqid)
                    qseq = String(row.qseq)
                    prefix, suffix = get(affix_dict, sseqid, ("", ""))
                    reference = String(Dict(DB)[sseqid])
                    trimmed, distance = trim_and_align_sequence(qseq, prefix, suffix, reference, stats,
                        min_quality=minquality, sseqid=sseqid)
                    return (trimmed, distance)
                end

                blast_clusters[:, :aln_qseq] = [r[1] for r in results]
                blast_clusters[:, :aln_mismatch] = [r[2] for r in results]

                # Core coverage filter
                db_lengths = Dict{String, Int}(name => length(seq) for (name, seq) in DB)
                blast_clusters[:, :db_length] = [get(db_lengths, String(row.sseqid), 0) for row in eachrow(blast_clusters)]
                blast_clusters[:, :corecov] = map(eachrow(blast_clusters)) do row
                    row.db_length > 0 ? length(row.aln_qseq) / row.db_length : 0.0
                end
                before_corecov = nrow(blast_clusters)
                filter!(x -> x.corecov >= min_corecov, blast_clusters)
                @info "Core coverage filter (>= $min_corecov): kept $(nrow(blast_clusters)) / $before_corecov"

                @info "Alignment stats: $(stats.total_attempts) attempts, $(stats.prefix_failures) prefix, $(stats.suffix_failures) suffix failures"
                if verbose
                    for (gene_name, messages) in stats.failed_genes
                        for msg in messages
                            @info "  Failed [$gene_name]: $msg"
                        end
                    end
                end
            end
        else
            blast_clusters[:, :aln_qseq] = blast_clusters[:, :qseq]
            blast_clusters[:, :aln_mismatch] = blast_clusters[:, :mismatch]
        end

        # Apply output filters
        min_fullcount = parsed_args["search"]["blast"]["minfullcount"]
        min_fullfrequency = parsed_args["search"]["blast"]["minfullfreq"]
        min_length = parsed_args["search"]["blast"]["length"]

        filter!(x -> x.full_count >= min_fullcount, blast_clusters)
        filter!(x -> x.full_frequency >= min_fullfrequency, blast_clusters)
        filter!(x -> length(x.qseq) >= min_length, blast_clusters)
        if !keep_failed
            filter!(x -> x.aln_mismatch >= 0, blast_clusters)
        end
        filter!(x -> x.aln_mismatch <= parsed_args["search"]["blast"]["maxdist"], blast_clusters)
        blast_clusters[:, :allele_name] = map(
            x -> (x.aln_mismatch == 0) ? x.sseqid : unique_name(x.sseqid, x.aln_qseq) * " Novel",
            eachrow(blast_clusters)
        )
        output = always_gz(parsed_args["search"]["blast"]["output"])
        @info "Saving discovered gene sequences in compressed $output file"
        CSV.write(output, blast_clusters, compress=true, delim='\t')
    end
end
