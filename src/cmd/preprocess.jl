function add_preprocess_args!(s)
        @add_arg_table! s["preprocess"] begin
            "demultiplex"
                help = "Demultiplex indexed plate libraries into a TSV with per-read metadata"
                action = :command
        end

        dm = s["preprocess"]["demultiplex"]
        @add_arg_table! dm begin
            "fastq"
                help = "Input FASTQ file with reads (single-end)"
                required = true
            "indices"
                help = "TSV with well indices/barcodes and case identifiers"
                required = true
            "output"
                help = "Output TSV (gz auto-enabled) with demultiplexed reads and metadata"
                required = true
        end

        add_arg_group!(dm, "Read filtering", "demux_filter")
        @add_arg_table! dm begin
            "-l", "--length"
                help = "Minimum read length to keep"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 200
            "--case-filter-regex"
                help = "Regex to keep only cases matching pattern (e.g., '[ACDERF]')"
                arg_type = String
        end

        add_arg_group!(dm, "Indexing and output", "demux_out")
        @add_arg_table! dm begin
            "-f", "--forwardarrayindex"
                help = "Name of the forward array index to use for demultiplexing (if present)"
                arg_type = String
                default = ""
            "-s", "--split"
                help = "Write per-case FASTQ files"
                action = :store_true
        end

    return s
end
