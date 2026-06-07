module Data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using DataFrames
    using Statistics
    using CSV
    using MD5
    using UnicodePlots

    export load_fasta, plotgenes, unique_name, sequence_hash, load_demultiplex
    export concatenate_columns, validate_types, get_ratio_threshold
    export barplot_if_available, histogram_if_available, heatmap_if_available, boxplot_if_available
    export round_floats!

    """
        round_floats!(df; digits=4) -> df

    Round every floating-point column of `df` to `digits` decimal places in place, so written
    tables carry readable values instead of full Float64 precision. Integer/string columns are
    left untouched; `Inf`/`NaN` pass through. Suitable for ratio/frequency columns (0–1 range);
    do not use on tables with tiny/huge floats like BLAST evalue/bitscore.
    """
    function round_floats!(df::DataFrame; digits::Int=4)
        for c in propertynames(df)
            col = df[!, c]
            eltype(col) <: AbstractFloat || continue
            df[!, c] = round.(col; digits=digits)
        end
        return df
    end

    function sequence_hash(seq; digits=4)
        "S" * lpad(string(parse(Int, bytes2hex(MD5.md5(seq))[(end-(digits-1)):end], base=16) % 10^digits), digits, '0')
    end

    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name, "_")))_$(sequence_hash(sequence, digits=digits))"
    end

    function write_fastq(path, records)
        FASTQ.Writer(open(path, "w")) do writer
            for record in records; write(writer, record); end
        end
    end

    function write_gz_fastq(path, records)
        FASTQ.Writer(GzipCompressorStream(open(path, "w"))) do writer
            for record in records; write(writer, record); end
        end
    end

    function write_fasta(path, records)
        FASTA.Writer(open(path, "w")) do writer
            for record in records; write(writer, record); end
        end
    end

    function write_gz_fasta(path, records)
        FASTA.Writer(GzipCompressorStream(open(path, "w"))) do writer
            for record in records; write(writer, record); end
        end
    end

    function load_demultiplex(path; limit=nothing)
        table = CSV.File(path, delim='\t', types=Dict(:case => String), limit=limit) |> DataFrame
        @assert all([name in names(table) for name in ["well", "case", "name", "genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        @info "Reading $(path)."
        return table
    end

    function validate_identifier(id)
        pattern = r"[A-Z].*\d.*\*(.*S.*\d.*|.*\d+$)"
        letters = collect(id)
        underscore = sum(map(x -> '_' == x, letters))
        star = sum(map(x -> '*' == x, letters))
        return occursin(pattern, id) && (underscore <= 1) && (star <= 1)
    end

    function validate_sequence(seq::String)
        allowed_chars = Set(['A', 'T', 'G', 'C'])
        return all(c -> c in allowed_chars, seq)
    end

    """
        register_fasta_entry!(name, seq, path, name_to_seq, seq_to_name)

    Record one `(allele name, sequence)` pair. Throws if `name` or `seq` was already seen with
    a conflicting partner, or if the exact pair is duplicated.
    """
    function register_fasta_entry!(
        name::AbstractString,
        seq::AbstractString,
        path::AbstractString,
        name_to_seq::Dict{String,String},
        seq_to_name::Dict{String,String},
    )
        name = String(name)
        seq = String(seq)
        if haskey(name_to_seq, name)
            other = name_to_seq[name]
            other != seq && throw(ErrorException(
                "FASTA \"$path\": allele name \"$name\" appears with different sequences"))
            throw(ErrorException(
                "FASTA \"$path\": duplicate entry for allele name \"$name\""))
        end
        if haskey(seq_to_name, seq)
            other = seq_to_name[seq]
            other != name && throw(ErrorException(
                "FASTA \"$path\": identical sequence appears under different names " *
                "\"$other\" and \"$name\""))
            throw(ErrorException(
                "FASTA \"$path\": duplicate sequence under allele name \"$name\""))
        end
        name_to_seq[name] = seq
        seq_to_name[seq] = name
        return nothing
    end

    function _fasta_reader(path::AbstractString)
        io = open(path, "r")
        stream = endswith(path, ".gz") ? GzipDecompressorStream(io) : io
        return FASTA.Reader(stream), stream, io
    end

    function _close_fasta_reader!(reader, stream, io, path::AbstractString)
        close(reader)
        endswith(path, ".gz") && close(stream)
        close(io)
    end

    """
        load_fasta(path; unique=true, validate_format=false)

    Load FASTA records as `(name, sequence)` tuples using the full header (`FASTA.description`).

    When `unique=true` (default), aborts on duplicate allele names, duplicate sequences, repeated
    rows, or one name mapping to multiple sequences / one sequence mapping to multiple names.

    When `validate_format=true`, also require non-empty names and sequences and IMGT-like
    identifiers (`validate_identifier`).
    """
    function load_fasta(path::AbstractString; unique::Bool=true, validate_format::Bool=false)
        records = Vector{Tuple{String,String}}()
        name_to_seq = Dict{String,String}()
        seq_to_name = Dict{String,String}()
        reader, stream, io = _fasta_reader(path)
        for record in reader
            name = String(FASTA.description(record))
            seq = String(FASTA.sequence(record))
            if validate_format
                length(seq) == 0 && throw(ErrorException("Empty sequence found: $name"))
                !validate_identifier(name) && throw(ErrorException("Invalid identifier found: $name"))
                length(name) == 0 && throw(ErrorException("Empty identifier found: $name"))
            end
            unique && register_fasta_entry!(name, seq, path, name_to_seq, seq_to_name)
            push!(records, (name, seq))
        end
        _close_fasta_reader!(reader, stream, io, path)
        return records
    end

    function process_fastq(callback::Function, path)
        if endswith(path, ".gz")
            stream = GzipDecompressorStream(open(path, "r"))
        else
            stream = open(path, "r")
        end
        reader = FASTQ.Reader(stream)
        for record in reader
            callback(record)
        end
        close(reader)
        endswith(path, ".gz") && close(stream)
    end

    """Unicode bar plot of labels → counts (no-op for empty input). Extra kwargs (e.g. color) pass through."""
    function barplot_if_available(labels, counts; kwargs...)
        isempty(labels) && return nothing
        println(UnicodePlots.barplot(labels, counts; kwargs...))
        return nothing
    end

    """Unicode boxplot of one or more named numeric series (no-op for empty input)."""
    function boxplot_if_available(labels, data; kwargs...)
        (isempty(labels) || isempty(data)) && return nothing
        println(UnicodePlots.boxplot(labels, data; kwargs...))
        return nothing
    end

    """Unicode histogram of a numeric collection (no-op for empty input)."""
    function histogram_if_available(values; kwargs...)
        isempty(values) && return nothing
        println(UnicodePlots.histogram(collect(values); kwargs...))
        return nothing
    end

    """Unicode heatmap of a matrix (e.g. a 4×L base-composition profile)."""
    function heatmap_if_available(matrix; kwargs...)
        println(UnicodePlots.heatmap(matrix; kwargs...))
        return nothing
    end

    function plotgenes(df)
        gene_counts = combine(groupby(df, :gene), nrow => :count)
        sort!(gene_counts, :count, rev=true)
        barplot_if_available(gene_counts.gene, gene_counts.count)
    end

    """
        concatenate_columns(row, col_names)

    Concatenate content from named columns into a single string.
    """
    function concatenate_columns(row, col_names)
        join(getproperty(row, Symbol(col_name)) for col_name in col_names)
    end

    """
        validate_types(types)

    Validate RSS type strings against allowed values.
    """
    function validate_types(types)
        allowed_types = ["heptamer", "spacer", "nonamer"]
        for type in types
            if !(type in allowed_types)
                error("Invalid type: $type. Allowed types are: heptamer, spacer, nonamer.")
            end
        end
        if isempty(types)
            error("At least one type must be specified.")
        end
    end

    """
        get_ratio_threshold(expect_dict, row; type="allele_ratio")

    Look up a per-allele or per-gene ratio threshold from expect_dict.
    Returns the threshold value (as Float64) if found, otherwise 0.0.
    """
    function get_ratio_threshold(expect_dict, row; type="allele_ratio")
        val = get(expect_dict, row.db_name, nothing)
        if val !== nothing
            @info "Applying $type >= $val for $(row.db_name)"
            return Float64(val)
        end
        val = get(expect_dict, row.gene, nothing)
        if val !== nothing
            @info "Applying $type >= $val for $(row.db_name)"
            return Float64(val)
        end
        return 0.0
    end
end
