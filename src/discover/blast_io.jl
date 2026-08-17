const columns = ["qseqid", "sseqid", "pident", "nident", "length", "mismatch", "gapopen", "qcovs", "qcovhsp", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore", "sstrand", "qseq"]

load_csv(path::String; delim::Char='\t') = CSV.File(path, delim=delim) |> DataFrame
select_columns(df::DataFrame, cols::Vector{Symbol}) = unique(select(df, cols))

"""Basename with BLAST/FASTA/TSV compound suffixes stripped (`foo.bar.fasta.gz` → `foo.bar`)."""
function file_stem(path::AbstractString)
    b = basename(String(path))
    m = match(r"^(.*?)(?:\.(?:tsv|fasta|fa|fna|blast|gz))+$", b)
    m !== nothing && return String(m.captures[1])
    stem, _ = splitext(b)
    return stem
end

function save_to_fasta(records::AbstractVector, output_file::String)
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
