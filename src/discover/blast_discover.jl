nogaps(s) = replace(s, '-' => "")

"Subject coverage of a BLAST hit: aligned span on the subject ÷ subject length, in (0, 1]."
subject_coverage(sstart::Integer, send::Integer, slen::Integer) = (abs(send - sstart) + 1) / slen

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

    cache_key = blast_cache_key(tsv_path; db_fasta=combined_db_fasta, args=args,
                                min_read_length=min_read_length)
    run_dir = blast_run_subdir(work_dir, cache_key)
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

    clusters = combine(groupby(blast_df, [:well, :case, :sseqid, :qseq, :mismatch]),
                        :qseqid => length => :full_count,
                        :scov => maximum => :scov)
    rename!(clusters, :mismatch => :blast_mismatch)
    @info "Clusters after grouping: $(nrow(clusters)) rows"
    verbose && CSV.write(joinpath(run_dir, "clusters.tsv"), clusters)

    before = nrow(clusters)
    mismatch_vals = Float64.(clusters.blast_mismatch)
    filter!(x -> x.blast_mismatch <= max_dist, clusters)
    stage_report("BLAST mismatch ≤ $max_dist (--max-blast-mismatch)", nrow(clusters), before;
                 values=mismatch_vals, histogram=true)
    verbose && CSV.write(joinpath(run_dir, "clusters-mismatch.tsv"), clusters)

    transform!(clusters, :sseqid => ByRow(x -> split(x, "*")[1]) => :gene)

    return sort(clusters, [:well, :case, :sseqid], rev=false)
end
