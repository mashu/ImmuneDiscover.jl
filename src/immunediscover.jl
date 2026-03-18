module immunediscover
    include("cli.jl")
    include("data.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("profile.jl")
    include("exact.jl")
    include("heptamer.jl")
    include("hsmm.jl")
    include("bwa.jl")
    include("blast.jl")
    include("keyedsets.jl")
    include("cooccurrence.jl")
    include("haplotype.jl")
    include("fasta.jl")
    include("merge.jl")
    include("table.jl")

    using .cli
    using .data
    using .demultiplex
    using .simulate
    using .profile
    using .exact
    using .heptamer
    using .HSMM
    using .bwa
    using .blast
    using .keyedsets
    using .cooccurrence
    using .haplotype
    using .fasta
    using .merge
    using .table

    using CSV
    using DataFrames
    using UnicodePlots
    using Glob
    using Statistics
    using DataStructures
    using FASTX
    using PrecompileTools

    export load_fasta, blast_discover

    # Backward compatibility: tests import immunediscover.association
    const association = cooccurrence

    """
        concatenate_columns(row, col_names)

    Helper function to concatenate content from columns
    """
    function concatenate_columns(row, col_names)
        concatenated_string = ""
        for col_name in col_names
            concatenated_string *= getproperty(row, Symbol(col_name))
        end
        concatenated_string
    end

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

    function get_ratio_threshold(expect_dict, row; type="allele_ratio")
        if row.db_name in keys(expect_dict)
            @info "Applying $type >= $(expect_dict[row.db_name]) for $(row.db_name)"
            return expect_dict[row.db_name]
        elseif row.gene in keys(expect_dict)
            @info "Applying $type >= $(expect_dict[row.gene]) for $(row.db_name)"
            return expect_dict[row.gene]
        else
            return 0
        end
    end

    """
        real_main(args=[])

    Main function for Julia
    """
    function real_main(args=[])
        parsed_args = parse_commandline(args)
        if parsed_args !== nothing
            if get(parsed_args, "%COMMAND%", "") == "demultiplex"
                demultiplex.handle_demultiplex(parsed_args, cli.always_gz)
            end

            if get(parsed_args, "%COMMAND%", "") == "search"
                subcmd = get(parsed_args["search"], "%COMMAND%", "")
                if subcmd == "heptamer"
                    heptamer.handle_heptamer(parsed_args, immunediscover, cli.always_gz)
                elseif subcmd == "exact"
                    exact.handle_exact(parsed_args, immunediscover, cli.always_gz)
                elseif subcmd == "hsmm"
                    HSMM.handle_hsmm(parsed_args)
                elseif subcmd == "blast"
                    blast.handle_blast(parsed_args, immunediscover, cli.always_gz)
                end
            end

            if get(parsed_args, "%COMMAND%", "") == "analyze"
                subcmd = get(parsed_args["analyze"], "%COMMAND%", "")
                if subcmd == "cooccurrence"
                    cooccurrence.handle_cooccurrence(parsed_args, cli.always_gz)
                elseif subcmd == "haplotype"
                    haplotype.handle_haplotype(parsed_args)
                elseif subcmd == "bwa"
                    bwa.handle_bwa(parsed_args, immunediscover, cli.always_gz)
                end
            end

            if get(parsed_args, "%COMMAND%", "") == "table"
                table.handle_table(parsed_args, immunediscover, cli.always_gz)
            end

            if get(parsed_args, "%COMMAND%", "") == "fasta"
                subcmd = get(parsed_args["fasta"], "%COMMAND%", "")
                if subcmd == "merge"
                    merge.handle_merge(parsed_args)
                elseif subcmd == "diff"
                    fasta.handle_fasta_diff(parsed_args, immunediscover)
                elseif subcmd == "hash"
                    fasta.handle_fasta_hash(parsed_args, immunediscover)
                end
            end
        end
    end

    """
        julia_main()::Cint

    Entry point for PackageCompiler. Top-level catch is acceptable at CLI boundary.
    """
    function julia_main()::Cint
        try
            real_main(ARGS)
        catch
            Base.invokelatest(Base.display_error, Base.catch_stack())
            return 1
        end
        return 0
    end

    @compile_workload begin
        parse_commandline(["--help"])
    end
end
