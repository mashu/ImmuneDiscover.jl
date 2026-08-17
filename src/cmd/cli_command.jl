# One singleton type per subcommand. Routing is by multiple dispatch on the concrete
# type (no Dict{String,Function}, no boxed closures), so each `run_command` method is
# its own specialization. The only run-time step is mapping the parsed (group, name)
# strings to the matching singleton once per invocation (command_for); everything after
# is statically dispatched. `run_command` is declared here and given methods by the top
# module, where the handlers (in the various submodules) are in scope.

abstract type Command end

struct PreprocessDemultiplex <: Command end
struct DiscoverBlast        <: Command end
struct DiscoverHsmm         <: Command end
struct DiscoverSelftest     <: Command end
struct SearchExact          <: Command end
struct SearchHeptamer       <: Command end
struct SearchBwa            <: Command end
struct AnalyzeCooccurrence  <: Command end
struct AnalyzeHaplotype     <: Command end
struct TableOuterjoin       <: Command end
struct TableLeftjoin        <: Command end
struct TableTransform       <: Command end
struct TableAggregate       <: Command end
struct TableUnique          <: Command end
struct TableSort            <: Command end
struct TableFilter          <: Command end
struct TableSelect          <: Command end
struct TableFasta           <: Command end
struct TableCollect         <: Command end
struct TableExclude         <: Command end
struct FastaMerge           <: Command end
struct FastaDiff            <: Command end
struct FastaHash            <: Command end

"(group, subcommand) path a Command is reached by on the command line."
cli_path(::PreprocessDemultiplex) = ("preprocess", "demultiplex")
cli_path(::DiscoverBlast)        = ("discover", "blast")
cli_path(::DiscoverHsmm)         = ("discover", "hsmm")
cli_path(::DiscoverSelftest)     = ("discover", "selftest")
cli_path(::SearchExact)          = ("search", "exact")
cli_path(::SearchHeptamer)       = ("search", "heptamer")
cli_path(::SearchBwa)            = ("search", "bwa")
cli_path(::AnalyzeCooccurrence)  = ("analyze", "cooccurrence")
cli_path(::AnalyzeHaplotype)     = ("analyze", "haplotype")
cli_path(::TableOuterjoin)       = ("table", "outerjoin")
cli_path(::TableLeftjoin)        = ("table", "leftjoin")
cli_path(::TableTransform)       = ("table", "transform")
cli_path(::TableAggregate)       = ("table", "aggregate")
cli_path(::TableUnique)          = ("table", "unique")
cli_path(::TableSort)            = ("table", "sort")
cli_path(::TableFilter)          = ("table", "filter")
cli_path(::TableSelect)          = ("table", "select")
cli_path(::TableFasta)           = ("table", "fasta")
cli_path(::TableCollect)         = ("table", "collect")
cli_path(::TableExclude)         = ("table", "exclude")
cli_path(::FastaMerge)           = ("fasta", "merge")
cli_path(::FastaDiff)            = ("fasta", "diff")
cli_path(::FastaHash)            = ("fasta", "hash")

const COMMANDS = (
    PreprocessDemultiplex(),
    DiscoverBlast(), DiscoverHsmm(), DiscoverSelftest(),
    SearchExact(), SearchHeptamer(), SearchBwa(),
    AnalyzeCooccurrence(), AnalyzeHaplotype(),
    TableOuterjoin(), TableLeftjoin(), TableTransform(), TableAggregate(),
    TableUnique(), TableSort(), TableFilter(), TableSelect(), TableFasta(),
    TableCollect(), TableExclude(),
    FastaMerge(), FastaDiff(), FastaHash(),
)

"""
    command_for(parsed_args) -> Command or absent

Resolve the parsed top-level group and its `%COMMAND%` subcommand to the matching
Command singleton (the single run-time mapping; dispatch is static thereafter).
"""
function command_for(parsed_args)
    group = String(get(parsed_args, "%COMMAND%", ""))
    isempty(group) && return absent
    return command_in_group(optional(get(parsed_args, group, nothing)), group)
end

command_in_group(::Absent, _) = absent
function command_in_group(block::Present, group)
    sub = String(get(block.value, "%COMMAND%", ""))
    for c in COMMANDS
        cli_path(c) == (group, sub) && return Present(c)
    end
    return absent
end

"Run a resolved command. Methods are defined by the top module (handlers live there)."
function run_command end

export Command, command_for, run_command
