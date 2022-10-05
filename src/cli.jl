module cli
    using ArgParse

    """
    Return path that ends with .gz or with appeneded extension
    """
    function always_gz(file_path)
        if endswith(file_path, ".gz")
            output = file_path
        else
            output = file_path*".gz"
        end
    end

    export parse_commandline
    """
    Handle command line
    """
    function parse_commandline()
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "0.1",
                            add_version = true,
                            epilog = "Mateusz Kaduk <mateusz.kaduk@ki.se>")
        @add_arg_table! s begin
            "demultiplex"
                help = "Demultiplex plate library"
                action = :command
            "trim"
                help = "Trim reads with a 5' and 3' profiles"
                action = :command
            end

        @add_arg_table! s["demultiplex"] begin
            "fastq"
                help = "FASTQ file with reads"
                required = true
            "indices"
                help = "TSV file with barcodes and cases"
                required = true
            "output"
                help = "TSV file to save demultiplex data"
                required = true
        end
        
        @add_arg_table! s["trim"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with aligned gene sequences to build trimming profile"
            required = true
        "output"
            help = "TSV file to save demultiplex data"
            required = true
        "-l", "--length"
            help = "Profile length"
            default = 20
            arg_type = Int
            range_tester = (x->x > 1)
        "-p", "--position"
            help = "An optional switch to save start,stop columns instead of trimmed_sequence"
            action = :store_true
        end

        return parse_args(s)
    end
end