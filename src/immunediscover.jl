module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("data.jl")
    include("profile.jl")
    include("trim.jl")

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    using .trim
    using CSV
    using DataFrames

    function real_main()
        parsed_args = parse_commandline()
        if parsed_args["%COMMAND%"] == "demultiplex"
            @info "Demultiplexing"
            table = immunediscover.demultiplex.demux(parsed_args["demultiplex"]["fastq"],
                                                     parsed_args["demultiplex"]["indices"])
            output = cli.always_gz(parsed_args["demultiplex"]["output"])
            CSV.write(output, table, compress=true, delim='\t')
            @info "Demultiplexed data saved in compressed $output file"
        end

        if parsed_args["%COMMAND%"] == "trim"
            @info "Trimming"
            table = CSV.File(parsed_args["trim"]["input"], delim='\t') |> DataFrame
            db = [r[2] for r in data.load_fasta(parsed_args["trim"]["fasta"])]
            len = parsed_args["trim"]["length"]

            prof_start, prof_stop = trim.trim_profiles(db, len)
            segments = trim.find_segments(table, prof_start, prof_stop)
            table[:,:prefix] .= [s[1] for s in segments]
            table[:,:center] .= [s[2] for s in segments]
            table[:,:suffix] .= [s[3] for s in segments]

            output = cli.always_gz(parsed_args["trim"]["output"])
            CSV.write(output, table, compress=true, delim='\t')
            @info "Trimmed data saved in compressed $output file"
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
