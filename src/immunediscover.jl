module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("data.jl")
    include("profile.jl")
    include("trim.jl")
    include("exact.jl")

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    using .trim
    using .exact
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
            pos = parsed_args["trim"]["position"]

            prof_start, prof_stop = trim.trim_profiles(db, len)
            ok, fragments, region = trim.find_segments(table.genomic_sequence, prof_start, prof_stop)

            subtable = table[ok,:]
            if pos
                start = [r[1] for r in region]
                stop = [r[2] for r in region]
                subtable[:,:start] = start
                subtable[:,:stop] = stop
            else
                subtable[:,:trimmed_sequence] = fragments
            end
            output = cli.always_gz(parsed_args["trim"]["output"])
            CSV.write(output, subtable, compress=true, delim='\t')
            @info "Trimmed data saved in compressed $output file"
        end

        if parsed_args["%COMMAND%"] == "exact"
            @info "Exact search"
            table = CSV.File(parsed_args["exact"]["tsv"], delim='\t') |> DataFrame
            db = data.load_fasta(parsed_args["exact"]["fasta"])
            counts_df = exact.exact_search(table, db)
            sort!(counts_df, [:case, :db_name])
            output = cli.always_gz(parsed_args["exact"]["output"])
            CSV.write(output, counts_df, compress=true, delim='\t')
            @info "Exact search data saved in compressed $output file"
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
