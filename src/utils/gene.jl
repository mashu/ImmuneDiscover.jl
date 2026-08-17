module Gene
    # Immunoglobulin / TCR segment identity. Owned here so search, simulate, haplotype, and
    # discover all dispatch on the same types instead of importing a search module for an enum.

    export GeneType, VGene, DGene, JGene
    export parse_gene_type, gene_string, gene_type_from_name, majority_gene

    abstract type GeneType end
    struct VGene <: GeneType end
    struct DGene <: GeneType end
    struct JGene <: GeneType end

    const GENE_TYPE_MAP = Dict{String,GeneType}("V" => VGene(), "J" => JGene(), "D" => DGene())

    """
        parse_gene_type(s) -> GeneType

    Convert a gene type string ("V", "J", "D") to a dispatch-ready type.
    """
    function parse_gene_type(s::AbstractString)
        gt = get(GENE_TYPE_MAP, String(s), nothing)
        gt === nothing && error("Invalid gene type: $s")
        return gt
    end

    gene_string(::VGene) = "V"
    gene_string(::JGene) = "J"
    gene_string(::DGene) = "D"

    const SEGMENT_CHAR_TO_TYPE = Dict('V' => VGene(), 'D' => DGene(), 'J' => JGene())
    const IG_TR_SEGMENT_PATTERN = r"(?:IG|TR)[A-Z]*([VDJ])(?:\d|\*)"

    """
        gene_type_from_name(name) -> GeneType or nothing

    Infer V/D/J segment type from an immunoglobulin or TCR gene name (e.g. `IGHV1-2`, `TRBV1-1`).
    Returns nothing for unrelated names (controls, housekeeping) that lack an IG/TR V/D/J segment.
    """
    function gene_type_from_name(name::AbstractString)
        m = match(IG_TR_SEGMENT_PATTERN, String(name))
        m === nothing && return nothing
        return get(SEGMENT_CHAR_TO_TYPE, only(m.captures[1]), nothing)
    end

    gene_tally(::VGene) = (1, 0, 0)
    gene_tally(::DGene) = (0, 1, 0)
    gene_tally(::JGene) = (0, 0, 1)

    """
        majority_gene(types) -> GeneType

    Segment type with the highest count. Ties prefer V, then D, then J.
    """
    function majority_gene(types::AbstractVector{<:GeneType})
        nv = nd = nj = 0
        for t in types
            a, b, c = gene_tally(t)
            nv += a
            nd += b
            nj += c
        end
        nv >= nd && nv >= nj && return VGene()
        nd >= nj && return DGene()
        return JGene()
    end
end
