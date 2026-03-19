module Exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds
    using FASTX
    using Statistics

    # ========================== GeneType dispatch hierarchy ==========================

    abstract type GeneType end
    struct VGene <: GeneType end
    struct DGene <: GeneType end
    struct JGene <: GeneType end

    """
        parse_gene_type(s) -> GeneType

    Convert a gene type string ("V", "J", "D") to a dispatch-ready type.
    """
    parse_gene_type(s::AbstractString) = s == "V" ? VGene() : s == "J" ? JGene() : s == "D" ? DGene() : error("Invalid gene type: $s")

    gene_string(::VGene) = "V"
    gene_string(::JGene) = "J"
    gene_string(::DGene) = "D"

    """
        gene_type_from_name(name) -> GeneType or nothing

    Infer gene type from an allele/gene name (e.g. "IGHV1-2" → VGene()).
    Returns nothing if no V/D/J is found.
    """
    function gene_type_from_name(name::AbstractString)
        occursin("V", name) && return VGene()
        occursin("D", name) && return DGene()
        occursin("J", name) && return JGene()
        return nothing
    end

    export GeneType, VGene, DGene, JGene, parse_gene_type, gene_type_from_name

    # ========================== Typed rows for border statistics ==========================

    const BORDER_ROW = NamedTuple{(:case, :gene, :matched_total, :accepted_total, :rejected_border, :rejected_ratio), Tuple{String, String, Int, Int, Int, Float64}}
    const BORDER_GENE_ROW = NamedTuple{(:gene, :mean_rejected_ratio, :matched_total, :rejected_border, :num_donors), Tuple{String, Float64, Int, Int, Int}}

    const LAST_BORDER_STATS = Ref(Vector{BORDER_ROW}())
    const LAST_BORDER_GENE_STATS = Ref(Vector{BORDER_GENE_ROW}())

    # ========================== RSS tuple filtering ==========================

    function filter_tuple_by_types(types_as_strings, mandatory_key; tuple_data)
        if mandatory_key ∉ ["prefix", "suffix"]
            error("mandatory_key must be either 'prefix' or 'suffix'")
        end
        full_tuple_keys_as_strings = ["prefix", "sequence", "heptamer", "spacer", "nonamer", "suffix"]
        selected_data = Dict("sequence" => tuple_data.sequence)
        if mandatory_key in keys(tuple_data)
            selected_data[mandatory_key] = getfield(tuple_data, Symbol(mandatory_key))
        end
        for key_str in full_tuple_keys_as_strings
            key_sym = Symbol(key_str)
            if (key_str in types_as_strings || key_str in ["sequence", mandatory_key]) && key_sym in keys(tuple_data)
                selected_data[key_str] = getfield(tuple_data, key_sym)
            end
        end
        return NamedTuple{Tuple(Symbol.(keys(selected_data)))}(values(selected_data))
    end

    # Filter flanks by gene type for RSS mode (dispatch replaces string ternary)
    filter_rss_flanks(rss, flanks, ::VGene) = filter_tuple_by_types(rss, "prefix", tuple_data=flanks)
    filter_rss_flanks(rss, flanks, ::JGene) = filter_tuple_by_types(rss, "suffix", tuple_data=flanks)
    filter_rss_flanks(rss, flanks, ::DGene) = flanks

    # ========================== Merge helpers ==========================

    function merge_min!(dst::Dict{String,Int}, src::Dict{String,Int})
        for (k, v) in src
            dst[k] = haskey(dst, k) ? min(dst[k], v) : v
        end
        return dst
    end

    function merge_counts!(dst::Dict{Tuple{String,String},Int}, src::Dict{Tuple{String,String},Int})
        for (k, v) in src
            dst[k] = get(dst, k, 0) + v
        end
        return dst
    end

    # ========================== Range validation ==========================

    function validate_range(start_pos, end_pos, seq_length)
        (start_pos < 1 || end_pos > seq_length || start_pos > end_pos) && error("Invalid range")
    end

    # ========================== extract_flanking — RSS mode (dispatch) ==========================

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::VGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
        return (prefix=prefix, sequence=seq,
                heptamer=gs[ep+1:min(ep+7,length(gs))],
                spacer=gs[ep+8:min(ep+30,length(gs))],
                nonamer=gs[ep+31:min(ep+39,length(gs))])
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::JGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (nonamer=gs[max(1,sp-39):max(1,sp-31)],
                spacer=gs[max(1,sp-30):max(1,sp-8)],
                heptamer=gs[max(1,sp-7):(sp-1)],
                suffix=gs[ep+1:min(ep+n,length(gs))],
                sequence=seq)
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ::Nothing)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (pre_nonamer=gs[max(1,sp-28):max(1,sp-20)],
                pre_spacer=gs[max(1,sp-19):max(1,sp-8)],
                pre_heptamer=gs[max(1,sp-7):(sp-1)],
                sequence=seq,
                post_heptamer=gs[ep+1:min(ep+7,length(gs))],
                post_spacer=gs[ep+8:min(ep+19,length(gs))],
                post_nonamer=gs[ep+20:min(ep+28,length(gs))])
    end

    # ========================== extract_flanking — Extension mode (dispatch) ==========================

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::VGene, n::Int, ext::Int)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
        return (prefix=prefix, sequence=seq, suffix=gs[ep+1:min(ep+ext,length(gs))])
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::JGene, n::Int, ext::Int)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (prefix=gs[max(1,sp-ext):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+n,length(gs))])
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, ::DGene, n::Int, ext::Int)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return (prefix=gs[max(1,sp-ext):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+ext,length(gs))])
    end

    # ========================== extract_flanking — Per-side extension (dispatch) ==========================

    function extract_flanking(gs::String, range::Tuple{Int,Int}, gt::GeneType, n::Int, ext::Int, prefE::Int, sufE::Int)
        sp, ep = range; validate_range(sp, ep, length(gs))
        seq = gs[sp:ep]
        return extract_perside(gs, sp, ep, seq, gt, n, prefE, sufE)
    end

    function extract_perside(gs, sp, ep, seq, ::VGene, n, prefE, sufE)
        prefix = sp > n ? gs[(sp-n):(sp-1)] : gs[1:(sp-1)]
        return (prefix=prefix, sequence=seq, suffix=gs[ep+1:min(ep+sufE,length(gs))])
    end
    function extract_perside(gs, sp, ep, seq, ::JGene, n, prefE, sufE)
        return (prefix=gs[max(1,sp-prefE):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+n,length(gs))])
    end
    function extract_perside(gs, sp, ep, seq, ::DGene, n, prefE, sufE)
        return (prefix=gs[max(1,sp-prefE):(sp-1)], sequence=seq, suffix=gs[ep+1:min(ep+sufE,length(gs))])
    end

    # ========================== String interface (used by tests and CLI) ==========================

    function extract_flanking(gs::String, range::Tuple{Int,Int}, gene_type::String, n::Int, extension::Union{Int,Nothing}=nothing)
        extract_flanking(gs, range, parse_gene_type(gene_type), n, extension)
    end

    function extract_flanking(gs::String, range::Tuple{Int,Int}, gene_type::String, n::Int,
                              extension::Union{Int,Nothing}, prefix_ext::Union{Int,Nothing}, suffix_ext::Union{Int,Nothing})
        gt = parse_gene_type(gene_type)
        extension === nothing && return extract_flanking(gs, range, gt, n, nothing)
        left_ext = prefix_ext === nothing ? extension : prefix_ext
        right_ext = suffix_ext === nothing ? extension : suffix_ext
        return extract_flanking(gs, range, gt, n, extension, left_ext, right_ext)
    end

    # ========================== extension_overlaps_border (dispatch) ==========================

    function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::VGene, ext::Int, border::Int)
        (border <= 0 || ext <= 0) && return false
        rbs = max(1, rl - border + 1)
        ee = min(ep + ext, rl)
        return (ep + 1 <= ee) && (ee >= rbs)
    end

    function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::JGene, ext::Int, border::Int)
        (border <= 0 || ext <= 0) && return false
        lbe = min(border, rl)
        es = max(1, sp - ext)
        return (es <= sp - 1) && (es <= lbe)
    end

    function extension_overlaps_border(sp::Int, ep::Int, rl::Int, ::DGene, ext::Int, border::Int)
        (border <= 0 || ext <= 0) && return false
        lbe = min(border, rl); rbs = max(1, rl - border + 1)
        les = max(1, sp - ext); ree = min(ep + ext, rl)
        return ((les <= sp-1) && (les <= lbe)) || ((ep+1 <= ree) && (ree >= rbs))
    end

    function extension_overlaps_border(sp::Int, ep::Int, rl::Int, gene_type::String, ext::Int, border::Int)
        extension_overlaps_border(sp, ep, rl, parse_gene_type(gene_type), ext, border)
    end

    # ========================== Calibration helpers (dispatch) ==========================

    function collect_allowed_lengths!(pre, suf, ::VGene, gb, ext, sp, ep, rl, border)
        rbs = max(1, rl - border + 1)
        push!(get!(suf, gb, Int[]), min(ext, max(0, rbs - 1 - ep)))
    end
    function collect_allowed_lengths!(pre, suf, ::JGene, gb, ext, sp, ep, rl, border)
        lbe = min(border, rl)
        push!(get!(pre, gb, Int[]), min(ext, max(0, sp - 1 - lbe)))
    end
    function collect_allowed_lengths!(pre, suf, ::DGene, gb, ext, sp, ep, rl, border)
        lbe = min(border, rl); rbs = max(1, rl - border + 1)
        push!(get!(pre, gb, Int[]), min(ext, max(0, sp - 1 - lbe)))
        push!(get!(suf, gb, Int[]), min(ext, max(0, rbs - 1 - ep)))
    end

    function safe_quantile(values::Vector{Int}, p::Float64, fallback::Int)
        isempty(values) && return fallback
        sorted = sort(values)
        idx = max(1, min(length(sorted), ceil(Int, (1 - p) * length(sorted))))
        return sorted[idx]
    end

    function assign_calibrated!(pgp, pgs, ::VGene, pre_v, suf_v, ext, ap)
        for (g, arr) in suf_v; pgs[g] = safe_quantile(arr, ap, ext); end
    end
    function assign_calibrated!(pgp, pgs, ::JGene, pre_v, suf_v, ext, ap)
        for (g, arr) in pre_v; pgp[g] = safe_quantile(arr, ap, ext); end
    end
    function assign_calibrated!(pgp, pgs, ::DGene, pre_v, suf_v, ext, ap)
        for (g, arr) in pre_v; pgp[g] = safe_quantile(arr, ap, ext); end
        for (g, arr) in suf_v; pgs[g] = safe_quantile(arr, ap, ext); end
    end

    # ========================== Border rejection (dispatch) ==========================

    function should_reject_border(::VGene, pgp, pgs, gb, ext, sp, ep, rl, border)
        rbs = max(1, rl - border + 1)
        return min(ep + get(pgs, gb, ext), rl) >= rbs
    end
    function should_reject_border(::JGene, pgp, pgs, gb, ext, sp, ep, rl, border)
        lbe = min(border, rl)
        return max(1, sp - get(pgp, gb, ext)) <= lbe
    end
    function should_reject_border(::DGene, pgp, pgs, gb, ext, sp, ep, rl, border)
        lbe = min(border, rl); rbs = max(1, rl - border + 1)
        return (max(1, sp - get(pgp, gb, ext)) <= lbe) || (min(ep + get(pgs, gb, ext), rl) >= rbs)
    end

    # ========================== Utility ==========================

    function get_ratio(expect_dict, row, ratio)
        row.db_name in keys(expect_dict) && (@info "Skipping allelic ratio filters for $(row.db_name) in case $(row.case)"; return 0)
        row.gene in keys(expect_dict) && (@info "Skipping gene ratio filters for $(row.db_name) in case $(row.case)"; return 0)
        return ratio
    end

    # ========================== exact_search ==========================

    function exact_search(table, query, gene; mincount=10, minratio=0.01, affix=13, rss=["heptamer", "spacer", "nonamer"], extension=nothing, N=10, raw=nothing, expect_dict=Dict{String,Float64}(), sequence_lookup=nothing, border::Int=0, adjust_per_gene_extension::Bool=false, adjust_percent::Float64=1.0)
        gt = parse_gene_type(gene)
        @info "Using mincount: $mincount, minratio: $minratio for genes: $gene"
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"

        per_gene_prefix = Dict{String,Int}()
        per_gene_suffix = Dict{String,Int}()
        if extension !== nothing && border > 0 && adjust_per_gene_extension
            @info "Calibrating per-gene extension targeting ≥ $(Int(round(adjust_percent*100)))% safe reads"
            p_cal = Progress(nrow(table))
            tmp = Folds.map(eachrow(table)) do row
                next!(p_cal)
                local pre = Dict{String,Vector{Int}}()
                local suf = Dict{String,Vector{Int}}()
                rl = length(row.genomic_sequence)
                @inbounds for (name, seq) in query
                    m = findfirst(seq, row.genomic_sequence)
                    if m !== nothing
                        gb = first(split(string(name), '*'))
                        collect_allowed_lengths!(pre, suf, gt, gb, extension, minimum(m), maximum(m), rl, border)
                    end
                end
                (pre=pre, suf=suf)
            end
            pre_vals = Dict{String,Vector{Int}}()
            suf_vals = Dict{String,Vector{Int}}()
            for r in tmp
                for (g, v) in r.pre; haskey(pre_vals, g) ? append!(pre_vals[g], v) : (pre_vals[g] = copy(v)); end
                for (g, v) in r.suf; haskey(suf_vals, g) ? append!(suf_vals[g], v) : (suf_vals[g] = copy(v)); end
            end
            assign_calibrated!(per_gene_prefix, per_gene_suffix, gt, pre_vals, suf_vals, extension, adjust_percent)
            @info "Per-gene extension calibrated for $(length(union(collect(keys(per_gene_prefix)), collect(keys(per_gene_suffix))))) genes"
        end

        p = Progress(nrow(table))
        result = Folds.map(eachrow(table)) do row
            next!(p)
            case = row.case; well = row.well
            matches = Vector{NamedTuple}()
            totals = Dict{Tuple{String,String},Int}()
            accepted = Dict{Tuple{String,String},Int}()
            @inbounds for (name, seq) in query
                m = findfirst(seq, row.genomic_sequence)
                if m !== nothing
                    sp = minimum(m); ep = maximum(m)
                    if extension !== nothing && border > 0
                        gb = first(split(string(name), '*'))
                        key = (string(case), gb)
                        totals[key] = get(totals, key, 0) + 1
                        if adjust_per_gene_extension
                            should_reject_border(gt, per_gene_prefix, per_gene_suffix, gb, extension, sp, ep, length(row.genomic_sequence), border) && continue
                        else
                            extension_overlaps_border(sp, ep, length(row.genomic_sequence), gt, extension, border) && continue
                        end
                        accepted[key] = get(accepted, key, 0) + 1
                    end
                    local flanks
                    if extension !== nothing && adjust_per_gene_extension && (border > 0)
                        gb = first(split(string(name), '*'))
                        flanks = extract_flanking(row.genomic_sequence, (sp, ep), gt, affix, extension, get(per_gene_prefix, gb, extension), get(per_gene_suffix, gb, extension))
                    else
                        flanks = extract_flanking(row.genomic_sequence, (sp, ep), gt, affix, extension)
                    end
                    filtered_flanks = extension !== nothing ? flanks : filter_rss_flanks(rss, flanks, gt)
                    meta = (well=string(well), case=string(case), db_name=string(name))
                    if extension !== nothing
                        push!(matches, Base.merge(meta, Base.merge(filtered_flanks, (prefix_len=length(filtered_flanks.prefix), suffix_len=length(filtered_flanks.suffix)))))
                    else
                        push!(matches, Base.merge(meta, filtered_flanks))
                    end
                end
            end
            (matches=matches, totals=totals, accepted=accepted)
        end

        valid_match_lists = [r.matches for r in result if !isempty(r.matches)]
        result_df = isempty(valid_match_lists) ? DataFrame() : DataFrame(reduce(vcat, valid_match_lists))

        if extension !== nothing && border > 0
            totals_all = Dict{Tuple{String,String},Int}()
            accepted_all = Dict{Tuple{String,String},Int}()
            for r in result
                !isempty(r.totals) && merge_counts!(totals_all, r.totals)
                !isempty(r.accepted) && merge_counts!(accepted_all, r.accepted)
            end
            if isempty(totals_all)
                LAST_BORDER_STATS[] = BORDER_ROW[]; LAST_BORDER_GENE_STATS[] = BORDER_GENE_ROW[]
            else
                stats_rows = BORDER_ROW[]
                tm = 0; tr = 0
                for (c, g) in union(collect(keys(totals_all)), collect(keys(accepted_all)))
                    mt = get(totals_all, (c, g), 0); a = get(accepted_all, (c, g), 0); r = mt - a
                    push!(stats_rows, (case=c, gene=g, matched_total=mt, accepted_total=a, rejected_border=r, rejected_ratio=(mt > 0 ? r/mt : 0.0)))
                    tm += mt; tr += r
                end
                @info "Border filter rejected $tr of $tm potential matches ($(tm>0 ? round(100*tr/tm;digits=2) : 0.0)%)"
                LAST_BORDER_STATS[] = stats_rows
                sr = Dict{String,Float64}(); sm = Dict{String,Int}(); srej = Dict{String,Int}(); nd = Dict{String,Int}()
                for row in stats_rows
                    g = row.gene; sr[g] = get(sr, g, 0.0) + row.rejected_ratio; sm[g] = get(sm, g, 0) + row.matched_total
                    srej[g] = get(srej, g, 0) + row.rejected_border; nd[g] = get(nd, g, 0) + 1
                end
                LAST_BORDER_GENE_STATS[] = [
                    (gene=g, mean_rejected_ratio=(nd[g]>0 ? sr[g]/nd[g] : 0.0), matched_total=sm[g], rejected_border=srej[g], num_donors=nd[g])
                    for g in keys(sr)]
            end
        else
            LAST_BORDER_STATS[] = BORDER_ROW[]; LAST_BORDER_GENE_STATS[] = BORDER_GENE_ROW[]
        end

        df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
        transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
        transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)

        if sequence_lookup !== nothing
            @info "Adding isin_db column based on reference FASTA"
            df[!, :isin_db] = map(row -> get(sequence_lookup, row.sequence, false) ? "" : "Novel", eachrow(df))
        end

        raw !== nothing && CSV.write(raw*".gz", result_df, delim='\t', compress=true)

        transform!(groupby(df, [:well, :case, :gene]), :full_count => (x->x./maximum(x)) => :full_ratio)
        transform!(groupby(df, [:well, :case, :gene]), :count => (x->x./maximum(x)) => :ratio)
        sort!(df, [:full_count, :count], rev=[true, true])
        udf = sort(unique(df),[:well, :case, :gene, :db_name, :sequence])
        filter!(row -> (row.full_count >= mincount) & (row.full_ratio >= get_ratio(expect_dict, row, minratio)), udf)
        filter!(row -> (row.count >= mincount) && (row.ratio >= get_ratio(expect_dict, row, minratio)), udf)

        priority_columns = ["well", "case", "gene", "db_name", "count", "full_count", "ratio", "full_ratio"]
        remaining_columns = setdiff(names(udf), priority_columns)
        udf[:, vcat(priority_columns, remaining_columns)]
        gdf = groupby(udf, [:well, :case, :gene, :db_name, :sequence])
        udf_indexed = transform(gdf, :well => (x -> 1:length(x)) => :flank_index)
        return filter(x->x.flank_index <= N, udf_indexed)
    end

    # ========================== Post-processing ==========================

    function transform_counts(group_df, name; count_col=:count)
        ref_row = filter(row -> startswith(row.db_name, name), group_df)
        ref_count = 1
        well, case = first(map(r->(r.well, r.case), eachrow(unique(group_df, [:well,:case]))))
        isempty(ref_row) ? (@warn "Reference name $name not found in well $well and case $case") : (@info "Applying name $name to well $well and case $case"; ref_count = ref_row.count)
        group_df[!, "$(count_col)_$(first(split(name,'*')))_ratio"] = group_df[:, count_col] ./ ref_count
        return group_df
    end

    function grouped_ratios(counts_df, refgene; count_col=:count)
        transformed = DataFrame[]
        for group in groupby(counts_df, [:well, :case])
            refgene != "" && (group = transform_counts(group, refgene, count_col=count_col))
            push!(transformed, DataFrame(group))
        end
        return reduce(vcat, transformed)
    end

    function build_sequence_lookup(ref_fasta_path::String)
        sequence_lookup = Dict{String, Bool}()
        @info "Building sequence lookup from reference FASTA: $ref_fasta_path"
        open(FASTA.Reader, ref_fasta_path) do reader
            for record in reader; sequence_lookup[string(FASTA.sequence(record))] = true; end
        end
        @info "Loaded $(length(sequence_lookup)) sequences from reference FASTA"
        return sequence_lookup
    end

    # ========================== CLI handler ==========================

    function handle_exact(parsed_args, immunediscover_module, always_gz)
        @info "Exact search"
        extension = parsed_args["search"]["exact"]["extension"]
        border = get(parsed_args["search"]["exact"], "border", 0)
        adjust_per_gene_extension = get(parsed_args["search"]["exact"], "adjust-per-gene-extension", false)
        adjust_percent = get(parsed_args["search"]["exact"], "adjust-percent", 1.0)
        limit = parsed_args["search"]["exact"]["limit"]
        refgenes = parsed_args["search"]["exact"]["refgene"]
        length(refgenes) > 0 && @info "Using reference genes $refgenes"

        table = immunediscover_module.load_demultiplex(parsed_args["search"]["exact"]["tsv"])
        limit > 0 && (@info "Limiting reads to $limit"; table = table[1:limit,:])
        db = immunediscover_module.load_fasta(parsed_args["search"]["exact"]["fasta"], validate=false)
        mincount = parsed_args["search"]["exact"]["mincount"]
        minratio = parsed_args["search"]["exact"]["minratio"]
        mincount < 5 && @warn "Decreasing mincount below 5 may lead to false positives"
        top = parsed_args["search"]["exact"]["top"]
        affix = parsed_args["search"]["exact"]["affix"]

        local rss
        if extension !== nothing
            @info "Using extension mode with length $extension"; rss = String[]
        else
            rss = split(parsed_args["search"]["exact"]["rss"], ',')
            immunediscover_module.validate_types(rss)
            @info "Extract RSS: $(join(rss,','))"
        end
        top != 1 && @info "Uncollapsed mode; at most $top full records returned."

        gene = parsed_args["search"]["exact"]["gene"]
        expect = parsed_args["search"]["exact"]["expect"]
        deletion = parsed_args["search"]["exact"]["deletion"]

        expect_df = DataFrame(name=[], ratio=[])
        expect !== nothing && (expect_df = CSV.read(expect, DataFrame, delim='\t'); @assert all([n in names(expect_df) for n in ["name","ratio"]]))
        nrow(expect_df) > 0 && @info "Using expect file with $(nrow(expect_df)) entries"
        expect_dict = Dict(zip(expect_df.name, expect_df.ratio))

        deletion_df = DataFrame(name=[], ratio=[])
        deletion !== nothing && (deletion_df = CSV.read(deletion, DataFrame, delim='\t'); @assert all([n in names(deletion_df) for n in ["name","ratio"]]))
        nrow(deletion_df) > 0 && @info "Using deletion file with $(nrow(deletion_df)) entries"
        deletion_dict = Dict(zip(deletion_df.name, deletion_df.ratio))

        raw = parsed_args["search"]["exact"]["raw"]
        locus = parsed_args["search"]["exact"]["locus"]
        ref_fasta = parsed_args["search"]["exact"]["ref-fasta"]
        sequence_lookup = ref_fasta !== nothing ? build_sequence_lookup(ref_fasta) : nothing

        counts_df = exact_search(table, db, gene, mincount=mincount, minratio=minratio,
            expect_dict=expect_dict, affix=affix, rss=rss, extension=extension, N=top,
            raw=raw, sequence_lookup=sequence_lookup, border=border,
            adjust_per_gene_extension=adjust_per_gene_extension, adjust_percent=adjust_percent)
        sort!(counts_df, [:case, :db_name])

        if !parsed_args["search"]["exact"]["noplot"]
            nrow(counts_df) > 0 ? immunediscover_module.plotgenes(counts_df) : @warn "No exact matches to plot"
        end

        @info "Excluding genes not starting with $locus for frequency calculation"

        transform!(groupby(counts_df, [:well, :case, :gene])) do group
            fg = filter(row -> startswith(row.db_name, locus), group)
            DataFrame(gene_count = fill(isempty(fg) ? 0 : sum(fg.count), nrow(group)))
        end
        transform!(groupby(counts_df, [:well, :case])) do group
            fg = filter(row -> startswith(row.db_name, locus), group)
            DataFrame(case_count = fill(isempty(fg) ? 0 : sum(fg.count), nrow(group)))
        end
        transform!(groupby(counts_df, [:gene])) do group
            fg = filter(row -> startswith(row.db_name, locus), group)
            DataFrame(cross_case_median_count = fill(isempty(fg) ? 0 : median(fg.count), nrow(group)))
        end
        transform!(groupby(counts_df, [:gene])) do group
            fg = filter(row -> startswith(row.db_name, locus), group)
            DataFrame(cross_case_median_gene_count = fill(isempty(fg) ? 0 : median(fg.gene_count), nrow(group)))
        end
        transform!(groupby(counts_df, [:db_name])) do group
            fg = filter(row -> startswith(row.db_name, locus), group)
            DataFrame(cross_case_median_allele_count = fill(isempty(fg) ? 0 : median(fg.count), nrow(group)))
        end

        counts_df[:,:allele_case_freq] = counts_df.count ./ counts_df.case_count
        counts_df[:,:gene_case_freq] = counts_df.gene_count ./ counts_df.case_count
        counts_df[:,:allele_to_cross_case_median_ratio] = counts_df.count ./ counts_df.cross_case_median_allele_count
        counts_df[:,:gene_to_cross_case_median_ratio] = counts_df.gene_count ./ counts_df.cross_case_median_gene_count
        transform!(groupby(counts_df, [:well, :case, :gene]), :count => (x->x./sum(x)) => :allele_freq)

        filter!(x -> x.allele_freq >= immunediscover_module.get_ratio_threshold(expect_dict, x, type="allele_freq"), counts_df)
        filter!(x -> x.gene_case_freq >= immunediscover_module.get_ratio_threshold(deletion_dict, x, type="gene_case_freq"), counts_df)

        mar = parsed_args["search"]["exact"]["min-allele-mratio"]
        mgr = parsed_args["search"]["exact"]["min-gene-mratio"]
        filter!(x -> x.allele_to_cross_case_median_ratio >= mar, counts_df)
        filter!(x -> x.gene_to_cross_case_median_ratio >= mgr, counts_df)

        output = always_gz(parsed_args["search"]["exact"]["output"])
        if length(refgenes) > 0
            for refgene in refgenes
                counts_df = grouped_ratios(counts_df, refgene, count_col=:count)
                transform!(groupby(counts_df, [:well, :case, :gene]), :count => sum => :ref_gene_count)
                counts_df = grouped_ratios(counts_df, refgene, count_col=:ref_gene_count)
            end
        end
        CSV.write(output, counts_df, compress=true, delim='\t')
        @info "Exact search data saved in compressed $output file"
    end

    export grouped_ratios, transform_counts, build_sequence_lookup, handle_exact
end
