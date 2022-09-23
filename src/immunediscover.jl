module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    using .cli
    using .demultiplex
    
    function real_main()
        parsed_args = parse_commandline()
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
