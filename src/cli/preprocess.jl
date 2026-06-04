function add_preprocess_args!(s)
        @add_arg_table! s["preprocess"] begin
            "demultiplex"
                help = "Demultiplex indexed plate libraries into a TSV with per-read metadata"
                action = :command
        end

        @add_arg_table! s["preprocess"]["demultiplex"] begin
            "fastq"
                help = "Input FASTQ file with reads (single-end)"
                required = true
            "indices"
                help = "TSV with well indices/barcodes and case identifiers"
                required = true
            "output"
                help = "Output TSV (gz auto-enabled) with demultiplexed reads and metadata"
                required = true
            "-l", "--length"
                help = "Minimum read length to keep"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 200
            "-s", "--split"
                help = "Write per-case FASTQ files"
                action = :store_true
            "-f", "--forwardarrayindex"
                help = "Name of the forward array index to use for demultiplexing (if present)"
                arg_type = String
                default = ""
            "--case-filter-regex"
                help = "Regex to keep only cases matching pattern (e.g., '[ACDERF]')"
                arg_type = String
        end

    return s
end
