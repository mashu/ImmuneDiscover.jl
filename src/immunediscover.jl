module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("data.jl")
    using .cli
    using .demultiplex
    using .simulate
    using .data
    using CSV

    function real_main()
        parsed_args = parse_commandline()
        if parsed_args["%COMMAND%"] == "demultiplex"
            @info "Demultiplexing"
            table = immunediscover.demultiplex.demux(parsed_args["demultiplex"]["fastq"], parsed_args["demultiplex"]["indices"])
            if endswith(parsed_args["demultiplex"]["output"],".gz")
                output = parsed_args["demultiplex"]["output"]
            else
                output = parsed_args["demultiplex"]["output"]*".gz"
            end
            CSV.write(output, table, compress=true, delim='\t')
            @info "Demultiplexed data saved in compressed $output file"
        end
        return 0    
    end    

    function julia_main()::Cint
        try
            real_main()
        catch
            Base.invokelatest(Base.display_error, Base.catch_stack())
            return 1
        end
        return 0
    end
end
