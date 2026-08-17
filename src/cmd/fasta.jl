function add_fasta_args!(s)
        @add_arg_table! s["fasta"] begin
            "merge"
                help = "Merge multiple FASTA files (aborts on duplicate names or sequences)"
                action = :command
            "diff"
                help = "Diff two FASTA files based on sequence identity but keep associated names"
                action = :command
            "hash"
                help = "Add hash based _S suffix to all allele names in the FASTA file"
                action = :command
        end

        mg = s["fasta"]["merge"]
        @add_arg_table! mg begin
        "output"
            help = "Output merged FASTA file"
            required = true
        "inputs"
            help = "Input FASTA files to merge (2 or more files)"
            nargs = '+'
            required = true
        end

        add_arg_group!(mg, "Naming", "merge_naming")
        @add_arg_table! mg begin
        "-c", "--cleanup"
            help = "Optional regex pattern to remove from sequence names (e.g., ' Novel')"
            arg_type = String
        "--add-source-prefix"
            help = "Add source filename as prefix to sequence names"
            action = :store_true
        end

        add_arg_group!(mg, "Merge behaviour", "merge_behaviour")
        @add_arg_table! mg begin
        "--no-sort"
            help = "Disable sorting sequences by name (default: sort enabled)"
            action = :store_true
        end

        @add_arg_table! s["fasta"]["diff"] begin
            "fasta"
                help = "FASTA files with sequences"
                nargs = '+'
                required = true
        end

        @add_arg_table! s["fasta"]["hash"] begin
        "fastain"
            help = "Input FASTA file path"
            required = true
        end

    return s
end
