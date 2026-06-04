# BLAST+ tabular output (outfmt 6 with custom fields): column layout, row parsing, and streaming reduction.
# Included into `Blast`; not a separate package module.

const BLAST_TABULAR_NAMES = (
    :qseqid, :sseqid, :pident, :nident, :length, :mismatch, :gapopen,
    :qcovs, :qcovhsp, :qstart, :qend, :sstart, :send, :qlen, :slen,
    :evalue, :bitscore, :sstrand, :qseq,
)
const columns = [string(n) for n in BLAST_TABULAR_NAMES]

const N_BLAST_TABULAR_FIELDS = length(BLAST_TABULAR_NAMES)

const BLAST_TABULAR_ELTYPE = (
    String, String, Float64, Int, Int, Int, Int, Int, Int, Int, Int, Int, Int, Int, Int,
    Float64, Float64, String, String,
)

blast_int(x::Missing) = 0
blast_int(x::Integer) = Int(x)
blast_int(x::AbstractFloat) = Int(round(x))
blast_int(x::AbstractString) = Int(round(parse(Float64, String(x))))
blast_int(x) = Int(round(Float64(x)))

blast_float(x::Missing) = 0.0
blast_float(x::Real) = Float64(x)
blast_float(x::AbstractString) = parse(Float64, String(x))
blast_float(x) = Float64(x)

blast_str(x::Missing) = ""
blast_str(x::AbstractString) = String(x)
blast_str(x) = String(x)

const BLAST_TABULAR_COERCE = (
    blast_str, blast_str, blast_float, blast_int, blast_int, blast_int, blast_int,
    blast_int, blast_int, blast_int, blast_int, blast_int, blast_int, blast_int, blast_int,
    blast_float, blast_float, blast_str, blast_str,
)

function blast_tabular_row_from_fields(fields::AbstractVector)
    length(fields) == N_BLAST_TABULAR_FIELDS ||
        error("expected $N_BLAST_TABULAR_FIELDS tab-separated BLAST fields, got $(length(fields))")
    NamedTuple{BLAST_TABULAR_NAMES}(
        ntuple(i -> BLAST_TABULAR_COERCE[i](fields[i]), N_BLAST_TABULAR_FIELDS),
    )
end

function blast_hit_preferred(a, b)
    (a.pident, a.qcovhsp, a.qcovs, a.bitscore) > (b.pident, b.qcovhsp, b.qcovs, b.bitscore)
end

function empty_blast_hits_dataframe()
    DataFrame([BLAST_TABULAR_NAMES[i] => BLAST_TABULAR_ELTYPE[i][] for i in 1:N_BLAST_TABULAR_FIELDS])
end

function best_hits_to_dataframe(best::Dict{String, <:NamedTuple})
    isempty(best) && return empty_blast_hits_dataframe()
    DataFrame(collect(values(best)); copycols=false)
end

"""Reposition `raw` at start and return either `raw` (plain TSV) or `BGZFStream(raw, \"r\")`. Caller must `close` only the returned stream."""
function blast_hits_decompressed_io(raw::IO)::IO
    mark(raw)
    g1 = read(raw, UInt8)
    g2 = read(raw, UInt8)
    if g1 != 0x1f || g2 != 0x8b
        reset(raw)
        return raw
    end
    reset(raw)
    return BGZFStream(raw, "r")
end

"""Open `path`, yield decompressed/tabular `io`; `close(io)` is always run after `f` (success or error)."""
function with_blast_hits_io(f::F, path::AbstractString) where {F}
    raw = open(path, "r")
    io = blast_hits_decompressed_io(raw)
    try
        return f(io)
    finally
        close(io)
    end
end

"""Parse each non-empty data line and pass the row to `visitor(row)` (typically a functor)."""
function scan_blast_tabular_lines(visitor::V, io::IO) where {V}
    for line in eachline(io)
        isempty(strip(line)) && continue
        fields = split(line, '\t'; keepempty=true)
        visitor(blast_tabular_row_from_fields(fields))
    end
    nothing
end

mutable struct BlastHitTableCollector
    rows::Vector{NamedTuple}
    BlastHitTableCollector() = new(NamedTuple[])
end

function (c::BlastHitTableCollector)(h::NamedTuple)
    push!(c.rows, h)
    c
end

mutable struct BlastBestHitPerQuery
    best::Dict{String, NamedTuple}
    n_raw::Int
    BlastBestHitPerQuery() = new(Dict{String, NamedTuple}(), 0)
end

function (sink::BlastBestHitPerQuery)(h::NamedTuple)
    sink.n_raw += 1
    qid = h.qseqid
    if !haskey(sink.best, qid) || blast_hit_preferred(h, sink.best[qid])
        sink.best[qid] = h
    end
    sink
end

"""
    read_blast_hits_table(path) -> DataFrame

Read all rows from a BLAST tabular file (plain TSV or BGZF `.gz` cache). Lines are parsed incrementally;
the full table is still held in the returned `DataFrame`.
"""
function read_blast_hits_table(path::AbstractString)::DataFrame
    collector = BlastHitTableCollector()
    with_blast_hits_io(path) do io
        scan_blast_tabular_lines(collector, io)
    end
    isempty(collector.rows) && return empty_blast_hits_dataframe()
    DataFrame(collector.rows; copycols=false)
end

"""
    stream_blast_best_hits(path) -> (DataFrame, n_raw::Int)

Single streaming pass over plain or **BGZF** BLAST tabular output: keep one row per `qseqid`
with the same winner rule as `combine(groupby(...), ...)`. Does not materialize all HSP rows in a `DataFrame`.
"""
function stream_blast_best_hits(path::AbstractString)::Tuple{DataFrame, Int}
    sink = BlastBestHitPerQuery()
    with_blast_hits_io(path) do io
        scan_blast_tabular_lines(sink, io)
    end
    best_hits_to_dataframe(sink.best), sink.n_raw
end
