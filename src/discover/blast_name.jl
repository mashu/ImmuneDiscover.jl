function seq_to_name_lookup(db_seqs)
    lookup = Dict{String, String}()
    for (name, seq) in db_seqs
        haskey(lookup, seq) || (lookup[seq] = name)
    end
    return lookup
end

"""
    name_candidate(core, sseqid, core_aln_mismatch, db_seqs, isin) -> String

Name a discovery candidate from its trimmed `core`. "Novel" means genuine sequence variation
in the (non-extended) gene relative to EVERY known allele, so the resolution order is:

  1. `core_aln_mismatch == 0`  → the best-hit reference `sseqid` (core equals that reference).
  2. core identical to ANY known sequence → that allele — independent of which allele won the
     BLAST best-hit and of `isin`. A core equal to a reference *is* that allele, never novel.
  3. `isin` only: core is a substring of a known allele → that allele (the read covers only
     part of the gene, with no internal variation — still known, not novel).
  4. otherwise → a hashed novel name (`unique_name`).

A core that strictly *contains* a known allele plus extra bases is deliberately left novel:
the extra bases are variation in the gene region (e.g. a junction insertion), not coverage.
`db_seqs` is the un-extended base reference as `(name, sequence)` pairs.
"""
function name_candidate(core::AbstractString, sseqid, core_aln_mismatch, db_seqs, isin::Bool)
    name_candidate(core, sseqid, core_aln_mismatch, seq_to_name_lookup(db_seqs), db_seqs, isin)
end

function name_candidate(core::AbstractString, sseqid, core_aln_mismatch,
                        exact_lookup::Dict{String, String},
                        db_seqs, isin::Bool)
    # Empty cores are trim failures, not known alleles. `occursin("", seq)` is true for every
    # reference, which previously named failed trims as the first DB allele.
    isempty(core) && return unique_name(sseqid, core)
    core_aln_mismatch == 0 && return String(sseqid)
    haskey(exact_lookup, core) && return exact_lookup[core]
    return name_from_substring(core, sseqid, db_seqs, Val(isin))
end

name_from_substring(core, sseqid, _, ::Val{false}) = unique_name(sseqid, core)
function name_from_substring(core, sseqid, db_seqs, ::Val{true})
    for (name, seq) in db_seqs
        isempty(seq) && continue
        occursin(core, seq) && return name
    end
    return unique_name(sseqid, core)
end
