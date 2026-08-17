module Gene
    # Immunoglobulin / TCR segment identity. Owned here so search, simulate, haplotype, and
    # discover all dispatch on the same types instead of importing a search module for an enum.

    export GeneType, VGene, DGene, JGene, Unsegmented
    export parse_gene_type, gene_string, gene_type_from_name, majority_gene, append_if_segmented!

    abstract type GeneType end
    struct VGene <: GeneType end
    struct DGene <: GeneType end
    struct JGene <: GeneType end
    struct Unsegmented <: GeneType end

    parse_gene_type(s::AbstractString) = parse_gene_type(Val{Symbol(s)}())
    parse_gene_type(::Val{:V}) = VGene()
    parse_gene_type(::Val{:D}) = DGene()
    parse_gene_type(::Val{:J}) = JGene()
    parse_gene_type(::Val{S}) where {S} = error("Invalid gene type: $S")

    gene_string(::VGene) = "V"
    gene_string(::JGene) = "J"
    gene_string(::DGene) = "D"
    gene_string(::Unsegmented) = ""

    segment_type(c::Char) = segment_type(Val(c))
    segment_type(::Val{'V'}) = VGene()
    segment_type(::Val{'D'}) = DGene()
    segment_type(::Val{'J'}) = JGene()
    segment_type(::Val{C}) where {C} = Unsegmented()

    const IG_TR_SEGMENT_PATTERN = r"(?:IG|TR)[A-Z]*([VDJ])(?:\d|\*)"

    """
        gene_type_from_name(name) -> GeneType

    Infer V/D/J segment type from an immunoglobulin or TCR gene name (e.g. `IGHV1-2`, `TRBV1-1`).
    Names without an IG/TR V/D/J segment (controls, housekeeping) are `Unsegmented`.
    """
    function gene_type_from_name(name::AbstractString)
        m = match(IG_TR_SEGMENT_PATTERN, String(name))
        m === nothing && return Unsegmented()
        return segment_type(only(m.captures[1]))
    end

    gene_tally(::VGene) = (1, 0, 0)
    gene_tally(::DGene) = (0, 1, 0)
    gene_tally(::JGene) = (0, 0, 1)
    gene_tally(::Unsegmented) = (0, 0, 0)

    append_if_segmented!(types, ::Unsegmented) = types
    append_if_segmented!(types, t::GeneType) = push!(types, t)

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
