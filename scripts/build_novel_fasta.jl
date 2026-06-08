#!/usr/bin/env julia
# Build Novel / Base / Complete FASTA sets from alleles.csv and SHORT references.
#
# Usage (from repo root):
#   julia --project=. scripts/build_novel_fasta.jl [output_dir]
#
# Reads:
#   data/alleles.csv
#   data/KI+1KGP-IGHV-SHORT.fasta
#   data/KI+1KGP-IGHD-SHORT.fasta
#
# Writes (default: data/):
#   KI+1KGP-IGHV-Novel.fasta, KI+1KGP-IGHD-Novel.fasta
#   KI+1KGP-IGHV-Base.fasta,   KI+1KGP-IGHD-Base.fasta
#   KI+1KGP-IGHV-Complete.fasta, KI+1KGP-IGHD-Complete.fasta

using CSV
using DataFrames
using FASTX

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = abspath(joinpath(SCRIPT_DIR, ".."))
const DEFAULT_DATA_DIR = joinpath(REPO_ROOT, "data")

const LOCUS_CONFIG = (
    (gene = "IGHV", short = "KI+1KGP-IGHV-SHORT.fasta", prefix = "IGHV"),
    (gene = "IGHD", short = "KI+1KGP-IGHD-SHORT.fasta", prefix = "IGHD"),
)

function clip_f_suffix(name::AbstractString)
    first(split(string(name), "_F"; limit=2))
end

function load_fasta_dict(path::AbstractString)
    records = Dict{String, String}()
    open(FASTA.Reader, path) do reader
        for record in reader
            name = string(strip(FASTA.identifier(record)))
            records[name] = string(FASTA.sequence(record))
        end
    end
    records
end

function write_fasta(path::AbstractString, records::Vector{Tuple{String, String}})
    open(FASTA.Writer, path) do writer
        for (name, seq) in records
            write(writer, FASTARecord(name, seq))
        end
    end
end

function novel_s_allele_names(table::DataFrame, gene_prefix::AbstractString)
    subtable = filter(row -> occursin(gene_prefix, row.db_name), table)
    subtable = filter(row -> row.KI_ImmuneDiscover > 0, subtable)
    unique_names = unique(string.(clip_f_suffix.(subtable.db_name)))
    filter(name -> occursin("_S", name), unique_names)
end

function sequence_for_name(
    name::AbstractString,
    short_db::Dict{String, String},
    allele_sequences::Dict{String, String},
)
    if haskey(short_db, name)
        return short_db[name]
    end
    if haskey(allele_sequences, name)
        @warn "Novel allele missing from SHORT FASTA; using alleles.csv sequence" name
        return allele_sequences[name]
    end
    error("No sequence found for novel allele: $name")
end

function build_locus_fasta(
    table::DataFrame,
    data_dir::AbstractString,
    gene::AbstractString,
    short_filename::AbstractString,
    gene_prefix::AbstractString,
)
    short_path = joinpath(data_dir, short_filename)
    short_db = load_fasta_dict(short_path)

    subtable = filter(row -> occursin(gene_prefix, row.db_name), table)
    novel_subtable = filter(row -> row.KI_ImmuneDiscover > 0, subtable)
    allele_sequences = Dict{String, String}()
    for row in eachrow(novel_subtable)
        key = clip_f_suffix(string(row.db_name))
        allele_sequences[key] = string(row.sequence)
    end

    novel_names = sort!(collect(String, novel_s_allele_names(table, gene_prefix)))
    novel_records = Tuple{String, String}[
        (name, sequence_for_name(name, short_db, allele_sequences))
        for name in novel_names
    ]

    novel_name_set = Set(first.(novel_records))
    base_records = Tuple{String, String}[
        (name, seq) for (name, seq) in sort(collect(short_db)) if name ∉ novel_name_set
    ]
    complete_records = vcat(novel_records, base_records)

    stem = "KI+1KGP-$(gene)"
    novel_path = joinpath(data_dir, "$(stem)-Novel.fasta")
    base_path = joinpath(data_dir, "$(stem)-Base.fasta")
    complete_path = joinpath(data_dir, "$(stem)-Complete.fasta")

    write_fasta(novel_path, novel_records)
    write_fasta(base_path, base_records)
    write_fasta(complete_path, complete_records)

    @info "Wrote $gene FASTA files" novel=length(novel_records) base=length(base_records) complete=length(complete_records)
end

function main()
    data_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : DEFAULT_DATA_DIR
    alleles_path = joinpath(data_dir, "alleles.csv")
    isfile(alleles_path) || error("Missing alleles table: $alleles_path")

    table = CSV.read(alleles_path, DataFrame, delim='\t')
    hasproperty(table, :KI_ImmuneDiscover) || error("alleles.csv must contain KI_ImmuneDiscover column")
    hasproperty(table, :db_name) || error("alleles.csv must contain db_name column")
    hasproperty(table, :sequence) || error("alleles.csv must contain sequence column")

    for config in LOCUS_CONFIG
        short_path = joinpath(data_dir, config.short)
        isfile(short_path) || error("Missing SHORT FASTA: $short_path")
        build_locus_fasta(table, data_dir, config.gene, config.short, config.prefix)
    end
end

main()
