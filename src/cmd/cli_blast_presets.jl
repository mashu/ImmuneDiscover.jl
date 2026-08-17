# `discover blast` tunable parameters: one source of truth.
# `BLAST_DEFAULTS` holds every tunable's global default; `discover.jl` reads them via
# `blast_default(key)` for its ArgParse `default=`, so the arg table and this table can
# never drift. Each gene preset lists ONLY the keys it changes from the default.

const BLAST_PARAM_GROUPS = [
    "Inputs / outputs"          => ["input", "fasta", "pseudo", "output", "full-output", "work-dir"],
    "Gene preset"               => ["gene", "show-presets"],
    "Extension & trimming"      => ["forward", "reverse", "minquality", "min-corecov"],
    "BLAST search"              => ["args", "max-blast-mismatch", "max-aln-mismatch", "edge", "subjectcov", "min-read-length"],
    "Cluster & output filters"  => ["length", "min-count", "min-fullcount", "min-allelic-ratio", "min-full-allelic-ratio", "min-peak-allelic-ratio", "min-reads-total", "min-recurrence", "max-homopolymer", "isin", "keep-failed"],
    "Run control"               => ["overwrite", "verbose"],
]

const BLAST_DEFAULTS = Dict{String,Any}(
    "forward"         => 20,
    "reverse"         => 20,
    "minquality"      => 0.75,
    "min-corecov"     => 0.6,
    "args"            => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50",
    "max-blast-mismatch" => 20,
    "max-aln-mismatch"   => 20,
    "edge"            => 0,
    "subjectcov"      => 0.1,
    "min-read-length" => 0,
    "min-count"       => 0,
    "min-fullcount"   => 5,
    "min-allelic-ratio"       => 0.0,
    "min-full-allelic-ratio"  => 0.1,
    "min-peak-allelic-ratio"  => 0.0,
    "min-reads-total" => 0,
    "length"          => 290,
    "min-recurrence"  => 0,
    "max-homopolymer" => 0,
    "work-dir"        => ".immunediscover",
)

"Global default for a `discover blast` parameter (single source for the ArgParse table)."
blast_default(key::AbstractString) = BLAST_DEFAULTS[key]

# Gene presets tuned on KI IGH self-tests (recovery of known-novel alleles; see selftest).
#   V: per-donor `--min-full-allelic-ratio` is off (0); the recall-safe cut is
#      `--min-peak-allelic-ratio 0.08` — a germline allele is a major allele (peak per-donor
#      allelic ratio ≥ 0.085) in at least one carrier, while PCR / sequencing artifacts never
#      are. 0.08 keeps every truth-novel allele (recall 1.0) while removing ~60% of false
#      novel calls vs the old 0.035. Leaving the global 0.1 per-donor floor on would drop
#      truth-novel alleles that peak above 0.08 in their best donor but sit below 0.1 there.
const BLAST_PRESETS = Dict(
    "V" => Dict{String,Any}(
        "min-full-allelic-ratio" => 0.0,
        "min-peak-allelic-ratio" => 0.08,
        "length"       => 283,
        "max-blast-mismatch" => 14,
        "max-aln-mismatch"   => 14,
        "minquality"   => 0.62,
        "min-corecov"  => 0.50,
    ),
    "D" => Dict{String,Any}(
        "forward"      => 40,
        "reverse"      => 40,
        "min-full-allelic-ratio" => 0.2,
        "length"       => 5,
        "min-fullcount" => 10,
        "edge"         => 10,
        "subjectcov"   => 0.25,
        "minquality"   => 0.5,
        "args"         => "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5",
    ),
    "J" => Dict{String,Any}(
        "forward"      => 12,
        "reverse"      => 12,
        "length"       => 10,
        "max-blast-mismatch" => 10,
        "max-aln-mismatch"   => 10,
        "min-fullcount" => 10,
        "args"         => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10",
    ),
)

"One parameter the gene preset wants to set; `applied` is false when a user override is kept."
struct PresetChange
    key::String
    current::Any
    preset::Any
    applied::Bool
end

"""
    preset_changes(block, gene) -> Vector{PresetChange}

Resolve what the `gene` preset would do to the parsed `block`. A preset value is applied
only when the current value still equals the global default (i.e. the user did not pass
it explicitly); otherwise the user override is kept. Sorted by key for stable logging.
"""
function preset_changes(block::AbstractDict, gene::AbstractString)
    preset = BLAST_PRESETS[gene]
    changes = PresetChange[]
    for key in sort!(collect(keys(preset)))
        haskey(block, key) || continue
        current = block[key]
        applied = current == blast_default(key)
        push!(changes, PresetChange(key, current, preset[key], applied))
    end
    return changes
end

"Single aligned block describing the preset outcome (replaces the per-key @info spam)."
function log_preset_changes(gene::AbstractString, changes::AbstractVector{PresetChange})
    printstyled("━━ $gene gene preset ", "━"^40, "\n"; color=:cyan, bold=true)
    applied = filter(c -> c.applied, changes)
    kept    = filter(c -> !c.applied, changes)
    kew = isempty(changes) ? 0 : maximum(length(c.key) for c in changes)
    if !isempty(applied)
        println("  applied (param was default → preset value)")
        for c in applied
            println("    ", rpad(c.key, kew), "  ", c.current, " → ", c.preset)
        end
    end
    if !isempty(kept)
        println("  kept (explicit user override; preset skipped)")
        for c in kept
            println("    ", rpad(c.key, kew), "  ", c.current, "  (preset ", c.preset, ")")
        end
    end
    return nothing
end

"Print every gene preset as the delta from the defaults (wired to `--show-presets`)."
function show_blast_presets()
    printstyled("BLAST gene presets (only keys that differ from the global default)\n";
                color=:cyan, bold=true)
    for gene in sort!(collect(keys(BLAST_PRESETS)))
        println("\n$gene gene:")
        preset = BLAST_PRESETS[gene]
        for key in sort!(collect(keys(preset)))
            println("  --$key  $(blast_default(key)) → $(preset[key])")
        end
    end
    return nothing
end

function get_blast_block(args)
    cmd = get(args, "%COMMAND%", "")
    cmd == "blast" && return Present(args["blast"])
    cmd == "discover" && get(args["discover"], "%COMMAND%", "") == "blast" &&
        return Present(args["discover"]["blast"])
    return absent
end

apply_blast_presets!(parsed_args, ::Absent) = parsed_args
function apply_blast_presets!(parsed_args, block::Present)
    b = block.value
    gene = b["gene"]
    haskey(BLAST_PRESETS, gene) || return parsed_args
    changes = preset_changes(b, gene)
    log_preset_changes(gene, changes)
    for c in changes
        c.applied && (b[c.key] = c.preset)
    end
    return parsed_args
end

"""
    apply_blast_presets!(parsed_args) -> parsed_args

Apply the selected gene preset in place: each preset key takes the preset value unless the
user passed it explicitly (detected as "current value ≠ default"). Logs one tidy block.
"""
apply_blast_presets!(parsed_args) = apply_blast_presets!(parsed_args, get_blast_block(parsed_args))
