module Fasta
    using CSV
    using DataFrames
    using Logging
    using ..Option: Absent, Present, absent, optional

    export extract_sequences_to_fasta, handle_fasta_diff, handle_fasta_hash

    include("fasta_extract.jl")
    include("fasta_handle.jl")
end
