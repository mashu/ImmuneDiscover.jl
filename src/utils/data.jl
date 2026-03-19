module Data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using DataFrames
    using Statistics
    using CSV
    using MD5
    using Requires

    const _barplot_fn = Ref{Union{Nothing, Function}}(nothing)
    function __init__()
        @require UnicodePlots = "b8865327-cd53-5732-bb35-84acbb429228" begin
            _barplot_fn[] = (x, y) -> UnicodePlots.barplot(x, y)
        end
    end

    export load_fasta, plotgenes, unique_name, sequence_hash, load_demultiplex
    export concatenate_columns, validate_types, get_ratio_threshold

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

    function load_fasta(path; validate=false)
        records = Vector{Tuple{String, String}}()
        seen_descriptions = Set{String}()
        seen_sequences = Set{String}()
        open(FASTA.Reader, path) do reader
            for record in reader
                desc = string(FASTA.description(record))
                seq = string(FASTA.sequence(record))
                if validate
                    length(seq) == 0 && throw(ErrorException("Empty sequence found: $(desc)"))
                    !validate_identifier(desc) && throw(ErrorException("Invalid identifier found: $(desc)"))
                    length(desc) == 0 && throw(ErrorException("Empty identifier found: $(desc)"))
                    desc in seen_descriptions && throw(ErrorException("Duplicate identifier found: $(desc)"))
                    seq in seen_sequences && throw(ErrorException("Duplicate sequence found: $(seq)"))
                end
                push!(seen_descriptions, desc)
                push!(seen_sequences, seq)
                push!(records, (desc, seq))
            end
        end
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

    """Print label/count as a simple text table when UnicodePlots is not available."""
    function print_counts_table(labels, counts)
        n = length(labels)
        n == 0 && return
        max_label = max(12, maximum(length(string(l)) for l in labels))
        max_count = maximum(length(string(c)) for c in counts)
        buf = IOBuffer()
        for i in 1:n
            println(buf, "  ", rpad(string(labels[i]), max_label), "  ", lpad(string(counts[i]), max_count))
        end
        println(String(take!(buf)))
    end

    """Barplot when UnicodePlots is available; otherwise print counts as a text table."""
    function barplot_if_available(x, y)
        f = _barplot_fn[]
        if f !== nothing
            println(f(x, y))
        else
            print_counts_table(x, y)
        end
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
    Returns the threshold value if found, otherwise 0.
    """
    function get_ratio_threshold(expect_dict, row; type="allele_ratio")
        val = get(expect_dict, row.db_name, nothing)
        if val !== nothing
            @info "Applying $type >= $val for $(row.db_name)"
            return val
        end
        val = get(expect_dict, row.gene, nothing)
        if val !== nothing
            @info "Applying $type >= $val for $(row.db_name)"
            return val
        end
        return 0
    end
end
