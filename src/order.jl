module alleleorder
    using CSV
    using DataFrames
    using DataStructures
    using FASTX
    using Logging
    using ProgressMeter
    using StatsBase
    using BurrowsWheelerAligner
    using CodecZlib
    using Base: @static
    using ..bwa

    """
        assign_reads_to_alleles(table::DataFrame, db::Vector{Tuple{String,String}}) -> DataFrame

    For each demultiplexed read, search for exact occurrence of any allele sequence.
    Return new DataFrame with an `allele_name` column (empty string if none matched)
    and `core_start`, `core_end` for matched substring positions (1-based, inclusive).
    """
    function assign_reads_to_alleles(table::DataFrame, db::Vector{Tuple{String,String}})
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        db_dict = Dict(db)
        allele_index = collect(keys(db_dict))
        p = Progress(nrow(table), desc="Assigning alleles by exact substring")
        allele_name = Vector{String}(undef, nrow(table))
        core_start = Vector{Int}(undef, nrow(table))
        core_end = Vector{Int}(undef, nrow(table))
        fill!(allele_name, "")
        fill!(core_start, 0)
        fill!(core_end, 0)
        for (i, row) in enumerate(eachrow(table))
            next!(p)
            seq = row.genomic_sequence
            matched = false
            for (name, aseq) in db
                m = findfirst(aseq, seq)
                if m !== nothing
                    allele_name[i] = String(name)
                    core_start[i] = minimum(m)
                    core_end[i] = maximum(m)
                    matched = true
                    break
                end
            end
            if !matched
                allele_name[i] = ""
                core_start[i] = 0
                core_end[i] = 0
            end
        end
        df = copy(table)
        df[:, :allele_name] = allele_name
        df[:, :core_start] = core_start
        df[:, :core_end] = core_end
        return df
    end

    """
        group_reads_by_allele(df::DataFrame; min_support::Int=5) -> Dict{String, Vector{NamedTuple}}

    Group full reads by `allele_name`, keeping those with at least `min_support` reads.
    Each group contains tuples (name, seq, core_start, core_end).
    """
    function group_reads_by_allele(df::DataFrame; min_support::Int=5, dedup::Bool=false, max_per_allele::Int=0)
        groups = Dict{String, Vector{NamedTuple}}()
        for g in groupby(df, :allele_name)
            a = String(first(g.allele_name))
            if isempty(a)
                continue
            end
            if nrow(g) < min_support
                continue
            end
            # Build read tuples with weights (counts) when dedup is enabled
            local reads
            if dedup
                # group by full sequence
                agg = combine(groupby(g, :genomic_sequence), nrow => :weight, first(:name) => :name, first(:core_start) => :core_start, first(:core_end) => :core_end)
                reads = [(name=String(r.name), seq=String(r.genomic_sequence), core_start=Int(r.core_start), core_end=Int(r.core_end), weight=Int(r.weight)) for r in eachrow(agg)]
            else
                reads = [(name=String(r.name), seq=String(r.genomic_sequence), core_start=Int(r.core_start), core_end=Int(r.core_end), weight=1) for r in eachrow(g)]
            end
            if max_per_allele > 0 && length(reads) > max_per_allele
                reads = reads[1:max_per_allele]
            end
            groups[a] = reads
        end
        return groups
    end

    """
        align_full_reads_to_genome(genome_paths::Vector{String}, reads_by_allele::Dict; chromosome::String, tag::String)
            -> Dict{String, Vector{NamedTuple}}

    Align full reads for each allele using BWA, keeping alignments on chromosome/tag.
    Returns positions per read with fields: ref, pos, mapq, nm, strand, core_exact::Bool.
    """
    function align_full_reads_to_genome(genome_paths::Vector{String}, reads_by_allele::Dict, db_core_by_allele::Dict{String,String}; chromosome::String, tag::String, external::Bool=false, threads::Int=4)
        if external
            return external_bwa_align(genome_paths, reads_by_allele, db_core_by_allele; threads=threads)
        end
        aligners = bwa.create_aligner(genome_paths)
        tag_regex = Regex(tag)
        results = Dict{String, Vector{NamedTuple}}()
        # Preload reference sequences for all genomes
        ref_seqs_dict = Dict{String, Dict{String, String}}()
        for genome_file in genome_paths
            ref_seqs_dict[genome_file] = bwa.load_reference_sequences(genome_file)
        end
        total_reads = sum(length(v) for v in values(reads_by_allele))
        p = Progress(total_reads, desc="Aligning full reads with BWA")
        for (allele, reads) in reads_by_allele
            placements = Vector{NamedTuple}()
            core_seq = get(db_core_by_allele, allele, "")
            for read in reads
                next!(p)
                record = FASTA.Record(read.name, read.seq)
                for (genome_file, aligner) in aligners
                    alns = BurrowsWheelerAligner.align(aligner, record)
                    if isempty(alns)
                        continue
                    end
                    # highest score alignments
                    max_score = maximum(map(x->x.score, alns))
                    alns = filter(x->x.score == max_score, alns)
                    descs = map(x->bwa.description(x, aligner), alns)
                    idx = map(x->occursin(chromosome, x) && occursin(tag_regex, x), descs)
                    if any(idx)
                        best = first(alns[idx])
                        ref = BurrowsWheelerAligner.refname(best, aligner)
                        strand = bwa.is_reverse_strand(best) ? "-" : "+"
                        # NM from packed field bits 10..31
                        nm = best.is_rev_is_alt_mapq_NM >> 10 & 0x003fffff
                        # Compute reference position for core start (0-based)
                        read_len = length(read.seq)
                        core_len = read.core_end - read.core_start + 1
                        if core_len <= 0
                            continue
                        end
                        core_ref_start0 = if strand == "+"
                            best.pos + (read.core_start - 1)
                        else
                            best.pos + (read_len - read.core_end)
                        end
                        # Extract reference sequence covering the core
                        ref_core_seq = bwa.extract_ref_subsequence(ref_seqs_dict[genome_file], String(ref), core_ref_start0, core_len)
                        if strand == "-"
                            ref_core_seq = bwa.reverse_complement_seq(ref_core_seq)
                        end
                        # Verify exact match of core to reference at aligned locus
                        core_exact = (!isempty(core_seq) && ref_core_seq == core_seq)
                        push!(placements, (read=read.name, ref=String(ref), pos=Int(core_ref_start0), nm=Int(nm), strand=String(strand), core_start=read.core_start, core_end=read.core_end, core_exact=core_exact, weight=read.weight))
                    end
                end
            end
            results[allele] = placements
        end
        return results
    end

    """
        external_bwa_align(genome_paths, reads_by_allele, db_core_by_allele; threads=4)

    Batch-align per-allele reads by writing temporary FASTA files and calling system bwa mem.
    We parse SAM and compute core_exact positions.
    """
    function external_bwa_align(genome_paths::Vector{String}, reads_by_allele::Dict, db_core_by_allele::Dict{String,String}; threads::Int=4)
        # Use first genome path for external bwa mem (assumes index available)
        genome = first(genome_paths)
        # Load reference sequences once for exact core check
        ref_seqs = bwa.load_reference_sequences(genome)
        results = Dict{String, Vector{NamedTuple}}()
        for (allele, reads) in reads_by_allele  
            # Create temp directory and files
            tmpdir = mktempdir()
            tmp_fa = joinpath(tmpdir, "reads.fa")
            open(tmp_fa, "w") do io
                for r in reads
                    println(io, ">", r.name)
                    println(io, r.seq)
                end
            end
            # Run bwa mem using Julia pipeline redirection
            sam = joinpath(tmpdir, "out.sam")
            open(sam, "w") do io
                run(pipeline(`bwa mem -t $(threads) $(genome) $(tmp_fa)`, stdout=io))
            end
            # Parse SAM and choose best placement per read (prefer core_exact, then MAPQ)
            placements = Vector{NamedTuple}()
            core_seq = get(db_core_by_allele, allele, "")
            # Accumulate candidates per read
            read_candidates = Dict{String, Vector{NamedTuple}}()
            open(sam, "r") do io
                for line in eachline(io)
                    startswith(line, "@") && continue
                    fields = split(line, '\t')
                    length(fields) >= 11 || continue
                    qname = fields[1]
                    flag = parse(Int, fields[2])
                    rname = fields[3]
                    pos1 = parse(Int, fields[4])
                    mapq = parse(Int, fields[5])
                    cigar = fields[6]
                    is_unmapped = (flag & 0x4) != 0
                    if is_unmapped
                        continue
                    end
                    is_rev = (flag & 0x10) != 0
                    # Lookup grouped read
                    matches = filter(r->r.name == qname, reads)
                    isempty(matches) && continue
                    read = first(matches)
                    core_len = read.core_end - read.core_start + 1
                    core_len <= 0 && continue
                    # Map read core interval to reference using CIGAR
                    core_start1, ok = map_read_pos_to_ref_start(cigar, pos1, read.core_start, read.core_end)
                    ok || continue
                    # Extract reference sequence at mapped core
                    ref_core_seq = bwa.extract_ref_subsequence(ref_seqs, String(rname), core_start1-1, core_len)
                    if is_rev
                        ref_core_seq = bwa.reverse_complement_seq(ref_core_seq)
                    end
                    core_exact = (!isempty(core_seq) && ref_core_seq == core_seq)
                    cand = (read=qname, ref=String(rname), pos=Int(core_start1-1), mapq=mapq, strand=is_rev ? "-" : "+", core_start=read.core_start, core_end=read.core_end, core_exact=core_exact, weight=read.weight)
                    push!(get!(read_candidates, qname, Vector{NamedTuple}()), cand)
                end
            end
            # Select best candidate per read: prefer core_exact, then highest mapq
            for (qname, cands) in read_candidates
                best = nothing
                # First choose those with core_exact
                exacts = filter(c->c.core_exact, cands)
                chosen = isempty(exacts) ? cands : exacts
                # pick max mapq
                best = findmax([(c.mapq, i) for (i,c) in enumerate(chosen)])[2]
                c = chosen[best]
                push!(placements, (read=c.read, ref=c.ref, pos=c.pos, nm=0, strand=c.strand, core_start=c.core_start, core_end=c.core_end, core_exact=c.core_exact, weight=c.weight))
            end
            rm(tmpdir; force=true, recursive=true)
            results[allele] = placements
        end
        return results
    end

    """
        map_read_pos_to_ref_start(cigar::String, pos1::Int, read_start::Int, read_end::Int) -> (ref_start1::Int, ok::Bool)

    Map a read interval [read_start, read_end] (1-based, inclusive) to the reference start (1-based)
    using the CIGAR string and leftmost reference position pos1. Returns ok=false if mapping fails.
    """
    function map_read_pos_to_ref_start(cigar::AbstractString, pos1::Int, read_start::Int, read_end::Int)
        # Parse CIGAR into operations
        ops = Vector{Tuple{Char,Int}}()
        num = ""
        for ch in cigar
            if '0' <= ch <= '9'
                num *= ch
            else
                isempty(num) && return 0, false
                push!(ops, (ch, parse(Int, num)))
                num = ""
            end
        end
        read_idx = 1
        ref_idx = pos1
        ref_core_start1 = -1
        # advance until we reach read_start
        for (op, len) in ops
            if op in ('M','X','=' )
                # both consume
                if read_idx + len - 1 < read_start
                    read_idx += len
                    ref_idx += len
                else
                    # read_start falls inside this block
                    delta = read_start - read_idx
                    ref_core_start1 = ref_idx + delta
                    break
                end
            elseif op in ('I','S')
                # consume read only
                if read_idx + len - 1 < read_start
                    read_idx += len
                else
                    # interval starts within an insertion/soft-clip - no reference mapping
                    return 0, false
                end
            elseif op in ('D','N')
                # consume ref only
                ref_idx += len
            elseif op in ('H','P')
                # consume neither
                continue
            else
                return 0, false
            end
        end
        if ref_core_start1 == -1
            # never reached start
            return 0, false
        end
        return ref_core_start1, true
    end

    """
        summarize_allele_positions(placements_by_allele::Dict) -> DataFrame

    Summarize per-allele genomic positions: count reads per (ref,pos) where core_exact is true.
    Also compute whether allele hits single or multiple loci.
    """
    function summarize_allele_positions(placements_by_allele::Dict)
        rows = Vector{NamedTuple}()
        for (allele, placements) in placements_by_allele
            if isempty(placements)
                continue
            end
            by_locus = Dict{Tuple{String,Int}, Int}()
            for p in placements
                if p.core_exact
                    key = (p.ref, p.pos)
                    w = hasproperty(p, :weight) ? p.weight :: Int : 1
                    by_locus[key] = get(by_locus, key, 0) + w
                end
            end
            total = sum(values(by_locus))
            for ((ref, pos), count) in by_locus
                push!(rows, (allele_name=String(allele), ref=String(ref), pos=Int(pos), count=Int(count), freq = total == 0 ? 0.0 : count/total))
            end
        end
        if isempty(rows)
            return DataFrame()
        end
        df = DataFrame(rows)
        sort!(df, [:ref, :pos, :allele_name])
        return df
    end

    """
        ascii_order(df::DataFrame) -> String

    Build a simple ASCII representation per reference of allele order along positions.
    """
    function ascii_order(df::DataFrame)
        if nrow(df) == 0
            return ""
        end
        out = IOBuffer()
        for ref_group in groupby(df, :ref)
            println(out, "== ", first(ref_group.ref), " ==")
            g = sort(ref_group, :pos)
            # Group by position and summarize overlapping alleles (bubble)
            parts = String[]
            for pos_group in groupby(g, :pos)
                alleles = sort(pos_group, :count, rev=true)
                label = join(["$(r.allele_name)[$(r.count)]" for r in eachrow(alleles)], "|")
                # Print 1-based position for readability
                push!(parts, "$(first(pos_group.pos)+1):{$(label)}")
            end
            println(out, join(parts, " -> "))
            println(out)
        end
        return String(take!(out))
    end

    """
        handle_order(parsed_args, immunediscover_module, always_gz)

    Coordinator for order command.
    """
    function handle_order(parsed_args, immunediscover_module, always_gz)
        @info "Allele genomic order inference"
        tsv = parsed_args["analyze"]["order"]["tsv"]
        fasta = parsed_args["analyze"]["order"]["fasta"]
        genome = parsed_args["analyze"]["order"]["genome"]
        output = parsed_args["analyze"]["order"]["output"]
        annotated = get(parsed_args["analyze"]["order"], "annotated", nothing)
        chromosome = parsed_args["analyze"]["order"]["chromosome"]
        tag = parsed_args["analyze"]["order"]["tag"]
        min_support = parsed_args["analyze"]["order"]["min-support"]
        dedup = parsed_args["analyze"]["order"]["dedup"]
        max_per_allele = parsed_args["analyze"]["order"]["max-per-allele"]
        external = parsed_args["analyze"]["order"]["external"]
        threads = parsed_args["analyze"]["order"]["threads"]
        ascii_out = get(parsed_args["analyze"]["order"], "ascii", nothing)

        demux = immunediscover_module.load_demultiplex(tsv)
        db = immunediscover_module.load_fasta(fasta, validate=false)

        assigned = assign_reads_to_alleles(demux, db)
        if annotated !== nothing
            CSV.write(always_gz(annotated), assigned, delim='\t', compress=true)
        end

        reads_by_allele = group_reads_by_allele(assigned, min_support=min_support, dedup=dedup, max_per_allele=max_per_allele)
        db_dict = Dict(db)
        placements_by_allele = align_full_reads_to_genome(genome, reads_by_allele, db_dict; chromosome=chromosome, tag=tag, external=external, threads=threads)
        summary_df = summarize_allele_positions(placements_by_allele)
        CSV.write(always_gz(output), summary_df, delim='\t', compress=true)
        @info "Order summary saved in $output"

        ascii_txt = ascii_order(summary_df)
        if ascii_out === nothing
            println(ascii_txt)
        else
            open(ascii_out, "w") do io
                write(io, ascii_txt)
            end
            @info "ASCII order visualization saved to $ascii_out"
        end
    end

    export handle_order
end


