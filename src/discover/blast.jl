module Blast
    using CSV
    using CodecZlib
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    using BioSequences
    using Folds
    using MD5
    using ProgressMeter: Progress, next!
    using Base.Threads: nthreads

    using ..Data: load_fasta as data_load_fasta, unique_name, histogram_if_available
    using ..SeqStats: gc_content, max_homopolymer
    using ..Filters: FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative,
                     add_group_ratio!, init_rejection_columns!, mark_rejected!, accepted, passes
    using ..Report: stage_report, section, cluster_profile_heatmap, params_report

    export blast_discover, save_to_fasta, accumulate_affixes, save_extended, handle_blast
    export resolve_work_dir, blast_hits_gz_path
    export consensus_prefix, consensus_suffix

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
        blastn(query_file, database, output_gz; args="")

    Run BLASTn: query sequences in `query_file` against BLAST DB at `database` (from `ensure_blast_db`).
    Streams tabular results to `output_gz` via gzip (`output_gz` must end in `.gz`; no uncompressed BLAST file is written).

    Thread count for the `blastn` process: `BLAST_NUM_THREADS` env if set (positive integer), else `Sys.CPU_THREADS`.
    (Julia's own thread pool is separate — start Julia with `-t auto` so trimming / `Folds` use multiple threads.)
    """
    function blastn_num_threads()
        v = strip(get(ENV, "BLAST_NUM_THREADS", ""))
        isempty(v) && return Sys.CPU_THREADS
        p = tryparse(Int, v)
        p === nothing && error("BLAST_NUM_THREADS must be a positive integer, got $(repr(v))")
        p < 1 && error("BLAST_NUM_THREADS must be >= 1, got $p")
        return p
    end

    # NCBI blastn expects single-dash long options (e.g. `-task`). GNU-style `--task` is rejected.
    function blastn_cli_token(t::AbstractString)::String
        s = String(t)
        if startswith(s, "--") && length(s) >= 3 && !startswith(s, "---")
            return string('-', SubString(s, 3))
        end
        return s
    end

    function build_blastn_cmd(query_file::String, database::String, outfmt::String, args::String)
        blast_num_threads = blastn_num_threads()
        cmd = `blastn -num_threads $blast_num_threads -query $query_file -db $database -out - -outfmt $outfmt`
        if !isempty(strip(args))
            # Julia 1.12+: `$(words...)` inside backticks concatenates into one argv; use `$words` for one arg each.
            extra = map(blastn_cli_token, split(args))
            cmd = `$cmd $extra`
        end
        return cmd
    end

    function blastn(query_file::String, database::String, output_gz::String; args::String="")
        endswith(output_gz, ".gz") || error("BLAST cache path must end with .gz, got $(repr(output_gz))")
        outfmt = "6 " * join(columns, " ")
        blast_num_threads = blastn_num_threads()
        cmd = build_blastn_cmd(query_file, database, outfmt, args)
        start_time = time()
        open(GzipCompressorStream, output_gz, "w") do gzio
            run(pipeline(cmd, stdout=gzio))
        end
        elapsed = time() - start_time
        @info "BLASTn completed in $(round(elapsed, digits=2)) seconds (blastn num_threads=$blast_num_threads, Julia nthreads=$(nthreads()), streamed to gzip)"
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

    """
        resolve_work_dir(path_str) -> String

    Absolute work directory for caches and intermediates. Relative paths are resolved against `pwd()`;
    empty `path_str` falls back to `joinpath(pwd(), \".immunediscover\")`.
    """
    function resolve_work_dir(path_str::AbstractString)::String
        s = strip(String(path_str))
        if isempty(s)
            return abspath(joinpath(pwd(), ".immunediscover"))
        end
        ex = expanduser(s)
        return isabspath(ex) ? abspath(ex) : abspath(joinpath(pwd(), ex))
    end

    """Stable subdirectory under `work_dir` for one demux input; used for BLAST cache and query
    FASTA. `salt` distinguishes runs whose query differs (e.g. a read-length filter)."""
    function blast_run_subdir(work_dir::AbstractString, input_tsv::AbstractString; salt::AbstractString="")::String
        tag = bytes2hex(md5(codeunits(abspath(input_tsv) * salt)))[1:16]
        return joinpath(work_dir, "blast", tag)
    end

    """Path to gzipped raw BLAST table for this input, matching `blast_discover` layout."""
    function blast_hits_gz_path(work_dir::AbstractString, input_tsv::AbstractString)::String
        d = blast_run_subdir(work_dir, input_tsv)
        return joinpath(d, "hits.blast.gz")
    end


    """Majority base at one flank column; ties broken lexicographically for stability."""
    function majority_base(counts::Dict{Char, Int})
        best_c = first(sort(collect(keys(counts)); by=c -> (-counts[c], c)))
        best_c => counts[best_c]
    end

    abstract type AffixSide end
    struct PrefixAffix <: AffixSide end
    struct SuffixAffix <: AffixSide end

    affix_column(a::AbstractString, col::Int, ::PrefixAffix) = a[length(a) - col + 1]
    affix_column(a::AbstractString, col::Int, ::SuffixAffix) = a[col]

    push_consensus_char!(out::Vector{Char}, c::Char, ::PrefixAffix) = pushfirst!(out, c)
    push_consensus_char!(out::Vector{Char}, c::Char, ::SuffixAffix) = push!(out, c)

    """
        consensus_affix(affixes, side; min_fraction) -> String

    Shared majority-vote extension from the gene boundary outward.
    `PrefixAffix`: 5' flanks, right-aligned, grow leftward.
    `SuffixAffix`: 3' flanks, left-aligned, grow rightward.
    """
    function consensus_affix(affixes::AbstractVector{<:AbstractString}, side::AffixSide;
                             min_fraction::Real=0.5)
        isempty(affixes) && return ""
        max_len = maximum(length, affixes)
        out = Char[]
        for col in 1:max_len
            counts = Dict{Char, Int}()
            n = 0
            for a in affixes
                length(a) >= col || continue
                c = affix_column(a, col, side)
                counts[c] = get(counts, c, 0) + 1
                n += 1
            end
            n == 0 && break
            best_c, best_n = majority_base(counts)
            best_n / n > min_fraction || break
            push_consensus_char!(out, best_c, side)
        end
        return String(out)
    end

    """
        consensus_prefix(prefixes; min_fraction) -> String

    Right-align upstream flanks (boundary-adjacent column first) and extend leftward
    while a strict majority of reads agree at each position.
    """
    consensus_prefix(prefixes::AbstractVector{<:AbstractString}; min_fraction::Real=0.5) =
        consensus_affix(prefixes, PrefixAffix(); min_fraction=min_fraction)

    """
        consensus_suffix(suffixes; min_fraction) -> String

    Left-align downstream flanks (boundary-adjacent column first) and extend rightward
    while a strict majority of reads agree at each position.
    """
    consensus_suffix(suffixes::AbstractVector{<:AbstractString}; min_fraction::Real=0.5) =
        consensus_affix(suffixes, SuffixAffix(); min_fraction=min_fraction)

    const AFFIX_ANCHOR_K = 12

    """Up to three k-mers per reference (5', middle, 3') for candidate filtering."""
    function gene_anchors(ref::AbstractString, k::Int)
        lr = length(ref)
        lr == 0 && return SubString{String}[]
        if lr <= k
            return [SubString(ref, 1, lr)]
        end
        mid = (lr - k) ÷ 2 + 1
        return unique([SubString(ref, 1, k), SubString(ref, mid, mid + k - 1), SubString(ref, lr - k + 1, lr)])
    end

    function build_affix_anchor_index(refs::AbstractVector{<:AbstractString}, k::Int)
        index = Dict{SubString{String}, Vector{Int}}()
        for (gi, ref) in enumerate(refs)
            for anchor in gene_anchors(ref, k)
                push!(get!(index, anchor, Int[]), gi)
            end
        end
        return index
    end

    function affix_match_positions(ref::AbstractString, gs::AbstractString)
        m = findfirst(ref, gs)
        m === nothing && return nothing
        return (minimum(m), maximum(m))
    end

    function affix_candidate_genes(gs::AbstractString, anchor_index::Dict{SubString{String}, Vector{Int}}, k::Int)
        gl = length(gs)
        gl == 0 && return Set{Int}()
        candidates = Set{Int}()
        if gl <= k
            for gi in get(anchor_index, SubString(gs, 1, gl), ())
                push!(candidates, gi)
            end
            return candidates
        end
        @inbounds for pos in 1:(gl - k + 1)
            anchor = SubString(gs, pos, pos + k - 1)
            for gi in get(anchor_index, anchor, ())
                push!(candidates, gi)
            end
        end
        return candidates
    end

    function affix_hits_for_read(gs::AbstractString, refs::AbstractVector{String},
                                 anchor_index::Dict{SubString{String}, Vector{Int}},
                                 forward_extension::Int, reverse_extension::Int)
        hits = Tuple{Int, String, String}[]
        gl = length(gs)
        for gi in affix_candidate_genes(gs, anchor_index, AFFIX_ANCHOR_K)
            pos = affix_match_positions(refs[gi], gs)
            pos === nothing && continue
            start_pos, end_pos = pos
            pre = start_pos > 1 ?
                  String(gs[max(1, start_pos - forward_extension):start_pos - 1]) : ""
            suf = end_pos < gl ?
                  String(gs[end_pos + 1:min(end_pos + reverse_extension, gl)]) : ""
            push!(hits, (gi, pre, suf))
        end
        return hits
    end

    function accumulate_affixes(db, demux_df; forward_extension=20, reverse_extension=20,
                                min_affix_fraction::Real=0.5)
        names = String[first(p) for p in db]
        refs = String[last(p) for p in db]
        n = length(refs)
        anchor_index = build_affix_anchor_index(refs, AFFIX_ANCHOR_K)
        prefixes = [String[] for _ in 1:n]
        suffixes = [String[] for _ in 1:n]

        prog = Progress(nrow(demux_df); desc="Collecting affixes")
        row_hits = Folds.map(eachrow(demux_df)) do row
            next!(prog)
            affix_hits_for_read(row.genomic_sequence, refs, anchor_index,
                                forward_extension, reverse_extension)
        end
        for hits in row_hits
            for (gi, pre, suf) in hits
                isempty(pre) || push!(prefixes[gi], pre)
                isempty(suf) || push!(suffixes[gi], suf)
            end
        end

        singleton = Vector{Tuple{String, String, String, String}}(undef, n)
        for gi in 1:n
            if isempty(prefixes[gi]) && isempty(suffixes[gi])
                singleton[gi] = (names[gi], refs[gi], "", "")
            else
                common_prefix = consensus_prefix(prefixes[gi]; min_fraction=min_affix_fraction)
                common_suffix = consensus_suffix(suffixes[gi]; min_fraction=min_affix_fraction)
                singleton[gi] = (names[gi], common_prefix * refs[gi] * common_suffix,
                                 common_prefix, common_suffix)
            end
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

    "Subject coverage of a BLAST hit: aligned span on the subject ÷ subject length, in (0, 1]."
    subject_coverage(sstart::Integer, send::Integer, slen::Integer) = (abs(send - sstart) + 1) / slen

    "Summarize cross-donor recurrence of accepted candidates (single-donor cores ≈ likely artifacts)."
    function report_recurrence(kept)
        nrow(kept) == 0 && return nothing
        uniq = unique(select(kept, [:aln_qseq, :n_donors]))
        n = nrow(uniq)
        single = count(==(1), uniq.n_donors)
        println("  $n distinct accepted candidate sequences; ",
                single, " seen in a single donor (", round(100 * single / n; digits=1),
                "% — more likely artifacts), ", n - single, " in two or more donors.")
        println("  donors per candidate (x = number of donors, bar height = candidates):")
        histogram_if_available(uniq.n_donors; nbins=20)
        return nothing
    end

    "Distance between two cores: Hamming when equal length, else Levenshtein."
    core_distance(a::AbstractString, b::AbstractString) =
        length(a) == length(b) ? sum(ca != cb for (ca, cb) in zip(a, b); init=0) :
                                 compute_edit_distance(String(a), String(b))

    """
        neighbor_stats(cores, reads; max_parents) -> (nn_dist, parent_ratio)

    For each distinct core, find its nearest *more-abundant* core (its likely "parent"): a
    small `nn_dist` with a large `parent_ratio` (= parent reads ÷ this core's reads) marks an
    error satellite of a dominant allele. Cores with no more-abundant neighbour get
    `nn_dist=-1`, `parent_ratio=1.0`. Only the `max_parents` most-abundant candidates are
    examined as potential parents (bounds cost; parents are abundant by definition).
    """
    function neighbor_stats(cores::AbstractVector{<:AbstractString}, reads::AbstractVector{<:Integer};
                            max_parents::Int=128)
        n = length(cores)
        nn = fill(-1, n)
        pr = ones(Float64, n)
        order = sortperm(reads, rev=true)            # most abundant first
        for ii in 1:n
            i = order[ii]
            best_d = typemax(Int)
            best_reads = 0
            considered = 0
            for jj in 1:(ii - 1)                      # candidates with ≥ reads
                j = order[jj]
                reads[j] > reads[i] || continue       # parent must be strictly more abundant
                considered += 1
                considered > max_parents && break
                d = core_distance(cores[i], cores[j])
                if d < best_d || (d == best_d && reads[j] > best_reads)
                    best_d, best_reads = d, reads[j]
                end
            end
            if best_reads > 0
                nn[i] = best_d
                pr[i] = best_reads / reads[i]
            end
        end
        return nn, pr
    end

    """
        add_neighbor_stats!(df; max_parents)

    Add `nn_dist` / `parent_ratio` columns (computed per gene over the distinct accepted cores;
    rejected rows keep the sentinels -1 / 1.0).
    """
    function add_neighbor_stats!(df::DataFrame; max_parents::Int=128)
        df[!, :nn_dist] = fill(-1, nrow(df))
        df[!, :parent_ratio] = ones(Float64, nrow(df))
        acc = df[df.reject_reason .== "", :]
        nrow(acc) == 0 && return df
        stats = Dict{String,Tuple{Int,Float64}}()       # core -> (nn_dist, parent_ratio)
        for gdf in groupby(acc, :gene)
            cores = String[]
            reads = Int[]
            seen = Set{String}()
            for r in eachrow(gdf)
                c = String(r.aln_qseq)
                c in seen && continue
                push!(seen, c); push!(cores, c); push!(reads, r.n_reads_total)
            end
            nd, prr = neighbor_stats(cores, reads; max_parents=max_parents)
            for k in eachindex(cores)
                stats[cores[k]] = (nd[k], prr[k])
            end
        end
        @inbounds for i in 1:nrow(df)
            df.reject_reason[i] == "" || continue
            s = get(stats, String(df.aln_qseq[i]), nothing)
            s === nothing && continue
            df.nn_dist[i], df.parent_ratio[i] = s
        end
        return df
    end

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
        # nt of read flanking the match on each side (symmetric: start-1 before, len-stop after)
        return five_prime - 1, length(read) - three_prime
    end

    """
        blast_discover(tsv_path, combined_db_fasta; kwargs...)

    Perform assignments and discovery of alleles based on BLAST results.
    """
    function blast_discover(tsv_path, combined_db_fasta; work_dir::AbstractString, max_dist=10, min_edge=10, min_scov=0.1, args="", verbose=false, overwrite=false, min_read_length::Int=0)
        min_ver = v"2.15.0"
        max_ver = v"2.17.0"
        is_valid, message = verify_blastn_version(min_ver, max_ver)
        @info message
        if !is_valid
            error("Please install BLAST version $min_ver - $max_ver")
        end

        # Fold the read-length threshold into the cache key so changing it doesn't reuse a
        # BLAST run built from a differently-filtered query.
        run_dir = blast_run_subdir(work_dir, tsv_path; salt = min_read_length > 0 ? "minlen=$min_read_length" : "")
        isdir(run_dir) || mkpath(run_dir)
        query_fasta = joinpath(run_dir, "query.fasta")
        blast_file = joinpath(run_dir, "hits.blast.gz")
        legacy_uncompressed = replace_extension(tsv_path, "blast")

        df = load_csv(tsv_path)
        @info "Read $(nrow(df)) rows from $tsv_path before BLAST assignment"
        unique_df = unique(df, :name)
        @info "Unique reads: $(nrow(unique_df))"
        if nrow(df) != nrow(unique_df)
            @error "Duplicated names found in input. Using unique ones but this is likely user error!"
        end
        df = unique_df

        # Drop short reads before BLAST (faster, less noise). Annotated into the cache key above.
        if min_read_length > 0
            before = nrow(df)
            filter!(r -> length(r.genomic_sequence) >= min_read_length, df)
            stage_report("read length ≥ $min_read_length nt (pre-BLAST)", nrow(df), before)
        end

        query_sequences = collect.(eachrow(select_columns(df, [:well, :case, :name, :genomic_sequence])))
        save_to_fasta(query_sequences, query_fasta)

        file_exists = isfile(blast_file)
        @info "BLAST cache check: file='$blast_file', exists=$file_exists, overwrite=$overwrite"
        if file_exists && !overwrite
            @info "BLASTn results already exist $blast_file. Skipping BLASTn."
        else
            @info "Running BLASTn. This may take a while."
            db_path = ensure_blast_db(combined_db_fasta)
            if isfile(legacy_uncompressed)
                rm(legacy_uncompressed)
                @info "Removed legacy uncompressed BLAST file beside input: $(basename(legacy_uncompressed))"
            end
            blastn(query_fasta, db_path, blast_file; args=args)
        end

        rm(query_fasta)

        blast_df = CSV.File(blast_file, delim='\t', header=columns) |> DataFrame
        @info "BLASTn results read from $blast_file: $(nrow(blast_df)) rows"
        if nrow(blast_df) == 0
            error("No BLASTn results found (wrong BLAST parameters?). Cannot proceed.")
        end

        # Keep only the single best hit per read. Pseudo genes (P-prefixed) are in the
        # database on purpose so they can WIN this best-hit step: a read whose best match is
        # a pseudo gene is then dropped below, excluding reads that look more like a pseudo
        # gene than a real allele. (Removing pseudo from the DB instead would let those reads
        # be assigned to their second-best real gene — a false positive.)
        blast_df = combine(groupby(blast_df, :qseqid), x -> first(sort(x, [:pident, :qcovhsp, :qcovs, :bitscore], rev=true)))
        leftjoin!(blast_df, df, on=:qseqid => :name)
        transform!(blast_df, [:qseq, :genomic_sequence] => ByRow(edge) => [:five_prime_edge, :three_prime_edge])

        @info "BLASTn results after best hits: $(nrow(blast_df)) rows"
        section("BLAST discovery — read & cluster filters")

        before = nrow(blast_df)
        filter!(x -> x.five_prime_edge >= min_edge && x.three_prime_edge >= min_edge, blast_df)
        stage_report("5'/3' edge ≥ $min_edge nt", nrow(blast_df), before)

        # Subject coverage from the sstart..send span (bounded ≤ 1), not BLAST's `length`
        # column which counts gap columns and can exceed slen on insertions.
        transform!(blast_df, [:sstart, :send, :slen] => ByRow(subject_coverage) => :scov)
        before = nrow(blast_df)
        scov_vals = copy(blast_df.scov)
        filter!(x -> x.scov > min_scov, blast_df)
        stage_report("subject coverage > $min_scov", nrow(blast_df), before; values=scov_vals, histogram=true, hist_closed=:right)

        transform!(blast_df, :qseq => ByRow(x -> replace(x, "-" => "")) => :qseq)
        before = nrow(blast_df)
        # Drop reads whose best hit was a pseudo gene (see best-hit note above).
        filter!(x -> !startswith(x.sseqid, "P"), blast_df)
        stage_report("drop pseudo-gene best hits", nrow(blast_df), before)
        verbose && CSV.write(joinpath(run_dir, "pseudo.tsv"), blast_df)

        read_name = Dict([(r.name, (r.well, r.case)) for r in eachrow(df)])
        blast_df[:, :well] = [read_name[x.qseqid][1] for x in eachrow(blast_df)]
        blast_df[:, :case] = [read_name[x.qseqid][2] for x in eachrow(blast_df)]

        clusters = combine(groupby(blast_df, [:well, :case, :sseqid, :qseq, :mismatch]), :qseqid => length => :full_count)
        @info "Clusters after grouping: $(nrow(clusters)) rows"
        verbose && CSV.write(joinpath(run_dir, "clusters.tsv"), clusters)

        before = nrow(clusters)
        mismatch_vals = Float64.(clusters.mismatch)
        filter!(x -> x.mismatch <= max_dist, clusters)
        stage_report("BLAST mismatch ≤ $max_dist", nrow(clusters), before; values=mismatch_vals, histogram=true)
        verbose && CSV.write(joinpath(run_dir, "clusters-mismatch.tsv"), clusters)

        transform!(clusters, :sseqid => ByRow(x -> split(x, "*")[1]) => :gene)
        add_group_ratio!(clusters, :full_count, [:well, :case, :gene], :full_ratio)

        return sort(clusters, [:well, :case, :sseqid], rev=false)
    end

    """
        handle_blast(parsed_args, immunediscover_module, always_gz)

    Handle the blast search pipeline. Extracted from real_main for modularity.
    """
    function handle_blast(parsed_args, immunediscover_module, always_gz)
        using_cli = immunediscover_module.Cli
        if get(parsed_args["discover"]["blast"], "show-presets", false)
            using_cli.show_blast_presets()
            return
        end
        parsed_args = using_cli.apply_blast_presets!(parsed_args)
        gene = parsed_args["discover"]["blast"]["gene"]
        haskey(using_cli.BLAST_PRESETS, gene) && @info "Applied $gene gene preset (overrides logged above)"
        # Grouped, ordered parameter display (final values after presets/overrides).
        params_report(parsed_args["discover"]["blast"], using_cli.BLAST_PARAM_GROUPS;
                      title="discover blast — parameters")
        @info "Discovery with BLAST assignments"

        fasta_path = parsed_args["discover"]["blast"]["fasta"]
        output_path = parsed_args["discover"]["blast"]["output"]
        work_dir = resolve_work_dir(parsed_args["discover"]["blast"]["work-dir"])
        isdir(work_dir) || mkpath(work_dir)
        @info "Work directory (caches, intermediates): $work_dir"
        file_stem = split(basename(fasta_path), '.')[1]
        affixes_path = joinpath(work_dir, file_stem * ".affixes")
        DB = immunediscover_module.load_fasta(fasta_path, validate=false)

        verbose = parsed_args["discover"]["blast"]["verbose"]
        overwrite = parsed_args["discover"]["blast"]["overwrite"]
        minquality = parsed_args["discover"]["blast"]["minquality"]
        forward_extension = parsed_args["discover"]["blast"]["forward"]
        reverse_extension = parsed_args["discover"]["blast"]["reverse"]
        keep_failed = parsed_args["discover"]["blast"]["keep-failed"]
        min_corecov = parsed_args["discover"]["blast"]["min-corecov"]
        isin = parsed_args["discover"]["blast"]["isin"]

        if (forward_extension < 7) && (forward_extension > 0)
            @warn "Forward extension $forward_extension is short and may lead to false positives"
        end
        if (reverse_extension < 7) && (reverse_extension > 0)
            @warn "Reverse extension $reverse_extension is short and may lead to false positives"
        end

        # Process pseudo genes if provided
        db_p = Vector{Tuple{String, String}}()
        pseudo = parsed_args["discover"]["blast"]["pseudo"]
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
                demux = load_csv(parsed_args["discover"]["blast"]["input"])
                @info "Extending gene sequences by $forward_extension forward and $reverse_extension reverse nucleotides" reads=nrow(demux) references=length(db_p) julia_threads=nthreads()
                # Extend the COMBINED db (real + pseudo) so the BLAST database matches the
                # no-extension branch. Pseudo decoys absent in reads become unextended
                # singletons but remain in the DB; using DB (real only) here silently dropped
                # them, disabling -p in the default extension mode.
                extended = accumulate_affixes(db_p, demux,
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
            parsed_args["discover"]["blast"]["input"],
            ext_fasta_path;
            work_dir=work_dir,
            max_dist=parsed_args["discover"]["blast"]["maxdist"],
            min_edge=parsed_args["discover"]["blast"]["edge"],
            min_scov=parsed_args["discover"]["blast"]["subjectcov"],
            args=parsed_args["discover"]["blast"]["args"],
            verbose=verbose,
            overwrite=overwrite,
            min_read_length=get(parsed_args["discover"]["blast"], "min-read-length", 0),
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
                section("BLAST discovery — trimming & core coverage")
                # corecov = trimmed core length (aln_qseq, the segment left after removing the
                # prefix/suffix affixes) ÷ length of the matched DB allele. Annotate, don't drop.
                mark_rejected!(blast_clusters, blast_clusters.corecov .< min_corecov,
                               "core coverage < $min_corecov (core÷DB len)", "core coverage")
                kept_corecov = count(isempty, blast_clusters.reject_reason)
                stage_report("core coverage ≥ $min_corecov (core÷DB len)", kept_corecov, before_corecov;
                             values=blast_clusters.corecov, histogram=true, hist_closed=:right)

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

        # Candidate-quality metrics (composition + cross-donor recurrence + support) — computed
        # before the output filters so they can also serve as optional filter criteria.
        #  - gc_content / max_homopolymer: composition of the trimmed core,
        #  - n_donors: distinct donors sharing the exact core (cross-donor recurrence),
        #  - n_reads_total: reads backing the core across the run.
        blast_clusters[:, :gc_content] = gc_content.(blast_clusters.aln_qseq)
        blast_clusters[:, :max_homopolymer] = max_homopolymer.(blast_clusters.aln_qseq)
        transform!(groupby(blast_clusters, :aln_qseq), :case => (x -> length(unique(x))) => :n_donors)
        transform!(groupby(blast_clusters, :aln_qseq), :full_count => sum => :n_reads_total)
        # Peak per-donor allelic ratio for the core: the highest fraction of its gene's reads
        # it reaches in any single donor. A germline allele is a major allele in ≥1 carrier;
        # artifacts stay minor everywhere — the most discriminative recall-safe separator.
        transform!(groupby(blast_clusters, :aln_qseq), :full_ratio => maximum => :max_full_ratio)

        # Apply output filters
        min_fullcount = parsed_args["discover"]["blast"]["minfullcount"]
        min_fullratio = parsed_args["discover"]["blast"]["minfullratio"]
        min_length = parsed_args["discover"]["blast"]["length"]
        min_recurrence = get(parsed_args["discover"]["blast"], "min-recurrence", 0)
        max_homop = get(parsed_args["discover"]["blast"], "max-homopolymer", 0)
        min_reads_total = get(parsed_args["discover"]["blast"], "min-reads-total", 0)

        criteria = FilterCriterion[
            MinThreshold(:full_count, min_fullcount, "min cluster reads (--minfullcount $min_fullcount)"),
            MinThreshold(:full_ratio, min_fullratio, "min allelic ratio (--minfullratio $min_fullratio)"),
            MinStringLength(:qseq, min_length, "min trimmed length (--length $min_length)"),
        ]
        if !keep_failed
            push!(criteria, NonNegative(:aln_mismatch, "trimming failed (aln_mismatch < 0)"))
        end
        push!(criteria, MaxThreshold(:aln_mismatch, Float64(parsed_args["discover"]["blast"]["maxdist"]),
                                     "max edit distance (--maxdist $(parsed_args["discover"]["blast"]["maxdist"]))"))
        # Optional quality-metric filters (off by default — 0 disables).
        min_reads_total > 0 && push!(criteria,
            MinThreshold(:n_reads_total, Float64(min_reads_total), "min total reads (--min-reads-total $min_reads_total)"))
        min_recurrence > 0 && push!(criteria,
            MinThreshold(:n_donors, Float64(min_recurrence), "min donor recurrence (--min-recurrence $min_recurrence)"))
        max_homop > 0 && push!(criteria,
            MaxThreshold(:max_homopolymer, Float64(max_homop), "max homopolymer (--max-homopolymer $max_homop)"))

        # Apply each output filter in turn, annotating (not dropping) so the full table keeps
        # every candidate, and report how many each filter removes — separately.
        section("BLAST discovery — output filters")
        init_rejection_columns!(blast_clusters)
        for criterion in criteria
            before = count(isempty, blast_clusters.reject_reason)
            fail = Bool[!passes(row, criterion) for row in eachrow(blast_clusters)]
            mark_rejected!(blast_clusters, fail, criterion.label, "output filter")
            stage_report(criterion.label, count(isempty, blast_clusters.reject_reason), before)
        end
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
            return unique_name(row.sseqid, row.aln_qseq)
        end
        blast_clusters[:, :allele_name] = map(allele_name_row, eachrow(blast_clusters))

        # Between-cluster separation (diagnostic): for each accepted core, the nearest
        # more-abundant core ("parent"). Small nn_dist + large parent_ratio = error satellite.
        add_neighbor_stats!(blast_clusters)

        # Two outputs: the filtered table (candidates that passed every stage) and a full table
        # holding every candidate plus reject_reason / reject_stage, for inspection and tuning.
        reason_cols = [:reject_reason, :reject_stage]
        output = always_gz(parsed_args["discover"]["blast"]["output"])
        full_arg = get(parsed_args["discover"]["blast"], "full-output", nothing)
        full_output = (full_arg === nothing || isempty(full_arg)) ?
            always_gz(replace(replace(output, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv") :
            always_gz(full_arg)

        section("BLAST discovery — summary")
        kept = accepted(blast_clusters)
        stage_report("accepted (passed all filters)", nrow(kept), nrow(blast_clusters))
        report_recurrence(kept)
        # Distance alone is not suspicious — most genuine novel alleles are 1 bp from a known
        # parent. What flags a likely error is being a small fraction of a much more abundant
        # neighbour, i.e. a large parent_ratio. Report both, but weight by parent_ratio.
        close_shadow = (kept.nn_dist .>= 0) .& (kept.nn_dist .<= 1) .& (kept.parent_ratio .>= 20)
        n_close = count((kept.nn_dist .>= 0) .& (kept.nn_dist .<= 1))
        n_shadow = count(close_shadow)
        n_close > 0 && @info "$n_close accepted candidate(s) are within 1 bp of a more-abundant core; $n_shadow of those carry < 5% of that neighbour's reads (parent_ratio ≥ 20) and are the more likely error satellites. nn_dist / parent_ratio are columns for your own threshold — not a default filter."
        # Quick look at how consistent the accepted candidates are (base composition over the
        # dominant length); no-op unless ≥2 accepted candidates share a length.
        nrow(kept) >= 2 && cluster_profile_heatmap(String.(kept.aln_qseq); title="accepted candidates")
        CSV.write(output, select(kept, Not(reason_cols)), compress=true, delim='\t')
        printstyled("  ✓ "; color=:green, bold=true); println("filtered  → $output  ($(nrow(kept)) rows)")
        CSV.write(full_output, blast_clusters, compress=true, delim='\t')
        printstyled("  ✓ "; color=:green, bold=true); println("full      → $full_output  ($(nrow(blast_clusters)) candidates + reject_reason)")
    end
end
