@testset "cli" begin
    @testset "version" begin
        @test Cli.software_version() == Cli.read_project_version()
        @test occursin(r"^\d+\.\d+\.\d+$", Cli.software_version())
        @test startswith(Cli.software_version_label(), Cli.software_version())
        @test occursin("git", Cli.software_version_label())
        @test Cli.cli_wants_version(["--version"])
        @test Cli.cli_wants_version(["-V"])
        @test !Cli.cli_wants_version(["search", "exact", "a.tsv", "b.fa", "c.tsv"])
        @test Cli.cli_is_help_or_version(["--help"])
        @test Cli.cli_is_help_or_version(["-h"])
        @test !Cli.cli_is_help_or_version(["search", "exact", "a.tsv", "b.fa", "c.tsv"])
    end

    @testset "command identity" begin
        paths = Set{Tuple{String,String}}()
        for cmd in Cli.COMMANDS
            group, sub = Cli.cli_path(cmd)
            @test !isempty(group)
            @test !isempty(sub)
            @test !occursin('-', group)
            @test !((group, sub) in paths)
            push!(paths, (group, sub))
            parsed = Dict{String,Any}(
                "%COMMAND%" => group,
                group => Dict{String,Any}("%COMMAND%" => sub),
            )
            resolved = Cli.command_for(parsed)
            @test resolved isa Present
            @test resolved.value === cmd
        end
        @test Cli.command_for(Dict{String,Any}()) === absent
        @test Cli.command_for(Dict{String,Any}(
            "%COMMAND%" => "search",
            "search" => Dict{String,Any}("%COMMAND%" => "missing"),
        )) === absent
    end

    @testset "help pages" begin
        settings = Cli.apply_parse_options!(Cli.build_argparse_settings(); exit_after_help=false)
        pages = Cli.cli_help_pages()
        @test [page.name for page in pages] == first.(Cli.argparse_named_nodes(settings))
        @test [page.args for page in pages] == Cli.help_invocation_args()

        help = mktemp() do path, io
            redirect_stdout(io) do
                Cli.parse_commandline(String["--help"]; exit_after_help=false)
            end
            flush(io)
            read(path, String)
        end
        @test help == Cli.help_text(settings, String["--help"])
        @test occursin("search", help)
        @test occursin("exact", help)
        @test occursin("--version", help)

        mktempdir() do help_dir
            Cli.write_cli_help_pages!(help_dir)
            names = Set(replace(f, ".txt" => "") for f in readdir(help_dir) if endswith(f, ".txt"))
            @test names == Set(page.name for page in pages)
            for page in pages
                page_path = joinpath(help_dir, page.name * ".txt")
                @test isfile(page_path)
                @test read(page_path, String) == Cli.help_text(settings, page.args)
            end
            for cmd in Cli.COMMANDS
                group, sub = Cli.cli_path(cmd)
                @test isfile(joinpath(help_dir, "$group.txt"))
                @test isfile(joinpath(help_dir, "$group-$sub.txt"))
            end
        end

        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                Cli.trace_cli_help_parse!()
            end
        end
    end

    @testset "run.sh help cache" begin
        Cli.ensure_help_pages!()
        run_sh = joinpath(@__DIR__, "..", "scripts", "run.sh")
        root = read(joinpath(Cli.CLI_HELP_DIR, "root.txt"), String)
        exact = read(joinpath(Cli.CLI_HELP_DIR, "search-exact.txt"), String)
        @test read(`$run_sh --help`, String) == root
        @test read(`$run_sh -h`, String) == root
        @test read(`$run_sh search exact --help`, String) == exact
        @test occursin("min-fullcount", exact)
        version = read(`$run_sh --version`, String)
        @test startswith(strip(version), Cli.software_version())
        @test occursin("git", version)
    end
end
