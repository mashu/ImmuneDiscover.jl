module cli
    using ArgParse
    
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
        
        return parse_args(s)
    end
end