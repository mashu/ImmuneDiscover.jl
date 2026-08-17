"""Return BLAST database path (no extension) for a FASTA. Build with makeblastdb if missing or stale."""
function ensure_blast_db(fasta_path::String)
    dir = dirname(fasta_path)
    db_path = joinpath(dir, file_stem(fasta_path))
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
    return joinpath(dirname(path), file_stem(path) * "." * ext)
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

"""
    blast_cache_key(input_tsv; db_fasta, args, min_read_length) -> String

Material hashed into `blast/<tag>/`. Must change when anything that affects raw BLASTn
output changes: query set, subject DB, blastn flags, or pre-BLAST read-length filter.
Gene preset / post-BLAST filters are intentionally excluded (they do not re-run blastn).
"""
function blast_cache_key(input_tsv::AbstractString;
                         db_fasta::AbstractString="",
                         args::AbstractString="",
                         min_read_length::Integer=0)
    parts = String[abspath(input_tsv)]
    !isempty(strip(String(db_fasta))) && push!(parts, "db=$(abspath(db_fasta))")
    !isempty(strip(String(args))) && push!(parts, "args=$(String(args))")
    min_read_length > 0 && push!(parts, "minlen=$min_read_length")
    return join(parts, "\0")
end

"""Stable subdirectory under `work_dir` for one BLAST cache key."""
function blast_run_subdir(work_dir::AbstractString, cache_key::AbstractString)::String
    tag = bytes2hex(md5(codeunits(String(cache_key))))[1:16]
    return joinpath(work_dir, "blast", tag)
end

"""Path to gzipped raw BLAST table, matching `blast_discover` layout."""
function blast_hits_gz_path(work_dir::AbstractString, input_tsv::AbstractString;
                            db_fasta::AbstractString="", args::AbstractString="",
                            min_read_length::Integer=0)::String
    key = blast_cache_key(input_tsv; db_fasta=db_fasta, args=args, min_read_length=min_read_length)
    return joinpath(blast_run_subdir(work_dir, key), "hits.blast.gz")
end

