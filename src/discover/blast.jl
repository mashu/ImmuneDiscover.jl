module Blast
    using CSV
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    using BioSequences
    using Folds

    using ..Data: load_fasta as data_load_fasta, unique_name
    using ..Filters: GermlineFilter, FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative, CustomFilter

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

    """Merge other into main (for combining thread-local stats after parallel map)."""
    function merge_stats!(main::AlignmentStats, other::AlignmentStats)
        main.total_attempts += other.total_attempts
        main.prefix_failures += other.prefix_failures
        main.suffix_failures += other.suffix_failures
        for (k, v) in other.successful_trims
            append!(main.successful_trims[k], v)
        end
        for (k, v) in other.failed_genes
            append!(main.failed_genes[k], v)
        end
        return main
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

    """Return BLAST database path (no extension) for a FASTA. Build with makeblastdb if missing or stale."""
    function ensure_blast_db(fasta_path::String)
        dir = dirname(fasta_path)
        base = first(split(basename(fasta_path), '.'))
        db_path = joinpath(dir, base)
        nin = db_path * ".nin"
        fasta_mtime = mtime(fasta_path)
        if isfile(nin) && mtime(nin) >= fasta_mtime
            @info "Using existing BLAST database $db_path"
            return db_path
        end
        @info "Building BLAST database from $fasta_path (enables multithreaded blastn)"
        run(`makeblastdb -in $fasta_path -dbtype nucl -parse_seqids -out $db_path`)
        return db_path
    end

    """
        blastn(query_file, database, output_file; ...)

    Run BLASTn: query sequences in `query_file` against `database`, write tabular results to `output_file`.
    With `use_db=true`, `database` is a BLAST DB path (from ensure_blast_db); multithreading is used.
    With `use_db=false`, `database` is a FASTA path and blastn uses -subject (single-threaded).
    """
    function blastn(query_file::String, database::String, output_file::String; args::String="", use_db::Bool=false)
        outfmt = "6 " * join(columns, " ")
        nthreads = Sys.CPU_THREADS
        if use_db
            cmd = `blastn -num_threads $nthreads -query $query_file -db $database -out $output_file -outfmt $outfmt`
        else
            cmd = `blastn -query $query_file -subject $database -out $output_file -outfmt $outfmt`
        end
        if !isempty(args)
            cmd = `$cmd $(split(args))`
        end
        elapsed = time_command(cmd)
        @info "BLASTn completed in $(round(elapsed, digits=2)) seconds" * (use_db ? " (num_threads=$nthreads)" : "")
    end

    """Map BLAST subject id to DB key. Novel alleles use base name (e.g. TRGV2*01_S2223 → TRGV2*01) for reference/affix lookup."""
    function sseqid_to_db_key(sseqid::AbstractString, db_keys::AbstractSet{<:AbstractString})
        s = String(sseqid)
        s in db_keys && return s
        m = match(r"^(.+)_S\d+$", s)
        if m !== nothing
            base = String(m.captures[1])
            base in db_keys && return base
        end
        return s
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

    # --- Affix trimming via dispatch (eliminates prefix/suffix code duplication) ---

    abstract type AffixSide end
    struct PrefixSide <: AffixSide end
    struct SuffixSide <: AffixSide end

    side_label(::PrefixSide) = "prefix"
    side_label(::SuffixSide) = "suffix"

    increment_failures!(stats, ::PrefixSide) = (stats.prefix_failures += 1)
    increment_failures!(stats, ::SuffixSide) = (stats.suffix_failures += 1)

    """
    For prefix: core is everything after the last aligned affix position.
    """
    function extract_core_range(aligned_affix, aligned_query, ::PrefixSide)
        boundary = findlast(x -> x != DNA_Gap, aligned_affix)
        boundary === nothing && return nothing, "No match found"
        start = boundary + 1
        start > length(aligned_query) && return nothing, "Query too short after trimming"
        return aligned_query[start:end], ""
    end

    """
    For suffix: core is everything before the first aligned affix position.
    """
    function extract_core_range(aligned_affix, aligned_query, ::SuffixSide)
        boundary = findfirst(x -> x != DNA_Gap, aligned_affix)
        boundary === nothing && return nothing, "No match found"
        boundary <= 1 && return nothing, "Starts too early"
        return aligned_query[1:boundary-1], ""
    end

    """
        trim_one_affix(affix, query, stats, scoremodel, side; min_quality, sseqid)

    Align `affix` to `query`, check quality, and extract the core sequence on the
    opposite side. Returns the trimmed query or nothing on failure.
    """
    function trim_one_affix(affix::LongDNA{4}, query::LongDNA{4}, stats, scoremodel, side::AffixSide;
                            min_quality=0.75, sseqid="")
        length(affix) == 0 && return query

        label = side_label(side)
        aln_result = safe_pairalign(affix, query, scoremodel)
        if aln_result === nothing
            push!(stats.failed_genes[sseqid], "$(titlecase(label)) alignment failed")
            increment_failures!(stats, side)
            return nothing
        end

        pairs = collect(alignment(aln_result))
        positions = findall(p -> first(p) != DNA_Gap, pairs)
        if isempty(positions)
            push!(stats.failed_genes[sseqid], "No $label content in alignment")
            increment_failures!(stats, side)
            return nothing
        end

        matches = sum(first(pairs[i]) == last(pairs[i]) for i in positions)
        quality = matches / length(positions)
        if quality < min_quality
            push!(stats.failed_genes[sseqid], "Poor $label alignment quality ($(round(quality * 100, digits=1))% match)")
            increment_failures!(stats, side)
            return nothing
        end

        aligned_affix = first.(pairs)
        aligned_query = last.(pairs)
        core, msg = extract_core_range(aligned_affix, aligned_query, side)
        if core === nothing
            push!(stats.failed_genes[sseqid], "$(titlecase(label)): $msg")
            increment_failures!(stats, side)
            return nothing
        end

        result = remove_gaps(LongDNA{4}(join(core)))
        if length(result) == 0
            push!(stats.failed_genes[sseqid], "Empty core after $label trimming")
            increment_failures!(stats, side)
            return nothing
        end
        return result
    end

    function trim_sequence(query::LongDNA{4}, prefix::LongDNA{4}, suffix::LongDNA{4}, stats,
        scoremodel::AffineGapScoreModel=AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1);
        min_quality=0.75, sseqid="")

        if length(query) == 0
            push!(stats.failed_genes[sseqid], "Empty query sequence")
            return nothing
        end
        stats.total_attempts += 1

        partial = trim_one_affix(prefix, query, stats, scoremodel, PrefixSide();
                                 min_quality=min_quality, sseqid=sseqid)
        partial === nothing && return nothing

        partial = trim_one_affix(suffix, partial, stats, scoremodel, SuffixSide();
                                 min_quality=min_quality, sseqid=sseqid)
        partial === nothing && return nothing

        if length(partial) < 10
            push!(stats.failed_genes[sseqid], "Core sequence too short ($(length(partial)) < 10 nt)")
            stats.prefix_failures += 1
            return nothing
        end

        push!(stats.successful_trims[sseqid], "Trimmed successfully")
        return partial
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
            db_path = ensure_blast_db(combined_db_fasta)
            blastn(query_fasta, db_path, blast_tsv, args=args, use_db=true)
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
        transform!(groupby(clusters, [:well, :case, :gene]), :full_count => (x -> x ./ maximum(x)) => :full_ratio)

        return sort(clusters, [:well, :case, :sseqid], rev=false)
    end

    """
        handle_blast(parsed_args, immunediscover_module, always_gz)

    Handle the blast search pipeline. Extracted from real_main for modularity.
    """
    function handle_blast(parsed_args, immunediscover_module, always_gz)
        using_cli = immunediscover_module.Cli
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
        output_path = parsed_args["search"]["blast"]["output"]
        work_dir = joinpath(dirname(abspath(output_path)), ".immunediscover")
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
        isin = parsed_args["search"]["blast"]["isin"]

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
            for (name, seq) in data_load_fasta(pseudo, validate=false)
                push!(db_p, ("P" * name, seq))
            end
        end
        for (name, seq) in data_load_fasta(fasta_path, validate=false)
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
                affix_dict = Dict{String, Tuple{String, String}}()
                for row in eachrow(CSV.File(affixes_path, delim='\t') |> DataFrame)
                    name = strip(ismissing(row.name) ? "" : String(row.name))
                    prefix = ismissing(row.prefix) ? "" : String(row.prefix)
                    suffix = ismissing(row.suffix) ? "" : String(row.suffix)
                    isempty(name) && continue
                    affix_dict[name] = (prefix, suffix)
                end

                stats = AlignmentStats()
                check_affix_quality_warning(forward_extension, minquality)
                check_affix_quality_warning(reverse_extension, minquality)

                # Normalize keys with strip() so FASTA headers with trailing space match BLAST sseqid (first token)
                db_keys = Set(strip(String(name)) for (name, seq) in DB)
                db_dict = Dict(strip(String(name)) => String(seq) for (name, seq) in DB)
                db_lengths = Dict(strip(String(name)) => length(seq) for (name, seq) in DB)

                # Parallel map with thread-local stats; merge into stats after (no concurrent mutation)
                results = Folds.map(eachrow(blast_clusters)) do row
                    local_stats = AlignmentStats()
                    sseqid = strip(String(row.sseqid))
                    base_key = sseqid_to_db_key(sseqid, db_keys)
                    qseq = String(row.qseq)
                    prefix, suffix = get(affix_dict, base_key, ("", ""))
                    reference = db_dict[base_key]
                    trimmed, distance = trim_and_align_sequence(qseq, prefix, suffix, reference, local_stats,
                        min_quality=minquality, sseqid=sseqid)
                    return (trimmed, distance, local_stats)
                end
                blast_clusters[:, :aln_qseq] = [r[1] for r in results]
                blast_clusters[:, :aln_mismatch] = [r[2] for r in results]
                for r in results
                    merge_stats!(stats, r[3])
                end

                # Core coverage filter (use base key for novel alleles with _S suffix)
                blast_clusters[:, :db_length] = [get(db_lengths, sseqid_to_db_key(strip(String(row.sseqid)), db_keys), 0) for row in eachrow(blast_clusters)]
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
        min_fullratio = parsed_args["search"]["blast"]["minfullratio"]
        min_length = parsed_args["search"]["blast"]["length"]

        criteria = FilterCriterion[
            MinThreshold(:full_count, min_fullcount, "Min cluster size"),
            MinThreshold(:full_ratio, min_fullratio, "Min allelic ratio"),
            MinStringLength(:qseq, min_length, "Min read length"),
        ]
        if !keep_failed
            push!(criteria, NonNegative(:aln_mismatch, "Exclude failed trimming"))
        end
        push!(criteria, MaxThreshold(:aln_mismatch, Float64(parsed_args["search"]["blast"]["maxdist"]), "Max distance"))
        GermlineFilter(criteria)(blast_clusters)
        # allele_name: exact match (aln_mismatch==0) -> sseqid; else if isin, substring of known allele -> that allele; else Novel
        db_seqs = [(strip(String(n)), String(s)) for (n, s) in DB]
        function allele_name_row(row)
            if row.aln_mismatch == 0
                return row.sseqid
            end
            if isin
                aln = String(row.aln_qseq)
                for (name, seq) in db_seqs
                    occursin(aln, seq) && return name
                end
            end
            return unique_name(row.sseqid, row.aln_qseq) * " Novel"
        end
        blast_clusters[:, :allele_name] = map(allele_name_row, eachrow(blast_clusters))
        output = always_gz(parsed_args["search"]["blast"]["output"])
        @info "Saving discovered gene sequences in compressed $output file"
        CSV.write(output, blast_clusters, compress=true, delim='\t')
    end
end
