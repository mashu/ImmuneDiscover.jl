module data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using DataFrames
    using Statistics
    using UnicodePlots
    using CSV
    using MD5

    export load_fasta, plotgenes, unique_name, sequence_hash, load_demultiplex

    function sequence_hash(seq; digits=4)
        "S"*lpad(string(parse(Int,bytes2hex(MD5.md5(seq))[(end-(digits-1)):end],base=16) % 10^digits), digits, '0')
    end

    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name,"_")))_$(sequence_hash(sequence,digits=digits))"
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
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        @info "Reading $(path)."
        return table
    end

    function validate_identifier(id)
        pattern = r"[A-Z].*\d.*\*(.*S.*\d.*|.*\d+$)"
        letters = collect(id)
        underscore = sum(map(x->'_' == x, letters))
        star = sum(map(x->'*' == x, letters))
        return occursin(pattern, id) && (underscore <= 1) && (star <= 1)
    end

    function validate_sequence(seq::String)
        allowed_chars = Set(['A', 'T', 'G', 'C'])
        return all(c -> c in allowed_chars, seq)
    end

    # FIX: Removed exit(1) calls — library code must throw, not terminate the process
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

    function plotgenes(df)
        gene_counts = combine(groupby(df, :gene), nrow => :count)
        sort!(gene_counts, :count, rev=true)
        barplot(gene_counts.gene, gene_counts.count)
    end
end
