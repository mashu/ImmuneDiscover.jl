function add_search_args!(s)
        @add_arg_table! s["search"] begin
            "exact"
                help = "Exact match to database alleles. Default TSV is slim (ids, counts, allelic ratios, flanks); pass --diagnostic for intermediate statistics."
                action = :command
            "heptamer"
                help = "Identify heptamer RSS positions and extend/trim V reads accordingly"
                action = :command
            "bwa"
                help = "Genome mapping QC with BWA to retain sequences mapped to the target chromosome"
                action = :command
        end

        bw = s["search"]["bwa"]
        @add_arg_table! bw begin
        "tsv"
            help = "TSV file with columns allele_name and seq"
            required = true
        "output"
            help = "TSV file to save filtered input"
            required = true
        "genome"
            help = "FASTA file with indexed genome"
            required = true
            nargs='+'
            arg_type = String
        end

        add_arg_group!(bw, "Target chromosome", "bwa_target")
        @add_arg_table! bw begin
        "-c", "--chromosome"
            help = "Chromosome string to filter by"
            default = "chromosome 14"
            arg_type = String
        "-t", "--tag"
            help = "Regex to filter valid descriptions of chromosomes"
            default = "(.*Primary Assembly.*)|(.*alternate locus.*)"
            arg_type = String
        end

        add_arg_group!(bw, "Input columns", "bwa_cols")
        @add_arg_table! bw begin
        "-n", "--colname"
            help = "Name of the column with allele names"
            default = "best_name"
            arg_type = String
        "-s", "--colseq"
            help = "List column names with sequences"
            default = ["prefix", "best_aln", "suffix"]
            nargs = '*'  # Accepts zero or more values
            arg_type = String
        end

        hp = s["search"]["heptamer"]
        @add_arg_table! hp begin
            "tsv"
                help = "TSV file with demultiplexed reads"
                required = true
            "fasta"
                help = "FASTA file with query alleles"
                required = true
            "output"
                help = "TSV file to save data with identified heptamers"
                required = true
            "summary"
                help = "TSV file to save summary collapsed alleles with statistics"
                required = true
        end

        add_arg_group!(hp, "Heptamer source", "hept_src")
        @add_arg_table! hp begin
            "-j", "--json"
                help = "JSON file with dictionary containing haptamers"
                default = "heptamers.json"
            "-c", "--chain"
                default = "IGHV"
                range_tester = (x->x ∈ CHAINS)
                arg_type = String
                help = "chain; must be one of " * join(CHAINS, ", ", " or ")
            "-d", "--maxdist"
                help = "A positive integer indicating maximum Hamming distance from any of heptamers in JSON file"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 1
        end

        add_arg_group!(hp, "Query trimming", "hept_trim")
        @add_arg_table! hp begin
            "-b", "--begin"
                help = "How much to trim from the 5' beginning of the query sequence"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 0
            "-e", "--end"
                help = "How much to trim from the 3' end of the query sequence"
                arg_type = Int
                range_tester = (x->x >= 0)
                default = 8
        end

        add_arg_group!(hp, "Summary filters", "hept_summary")
        @add_arg_table! hp begin
            "-m", "--mincount"
                help = "Minimum count allowed in summary"
                default = 1
                arg_type = Int
                range_tester = (x->x >= 1)
            "-r", "--ratio"
                help = "Lowest allowed ratio between counts of full allele sequence and trimmed allele sequence"
                default = 0.25
                arg_type = Float64
                range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        ex = s["search"]["exact"]

        @add_arg_table! ex begin
        "tsv"
            help = "TSV file with demultiplexed data"
            required = true
        "fasta"
            help = "FASTA file with query alleles"
            required = true
        "output"
            help = "TSV file for filtered results (a .full.tsv.gz sibling lists every candidate + reject_reason)"
            required = true
        end

        add_arg_group!(ex, "Gene selection", "exact_gene")
        @add_arg_table! ex begin
        "-g", "--gene"
            default = "V"
            range_tester = (x->x ∈ GENES)
            arg_type = String
            help = "gene; must be one of " * join(GENES, ", ", " or ")
        "--locus"
            help = "Optional db_name prefix for locus-scoped frequency columns (e.g. IGHV, TRBV, IG). Empty = all alleles included. Non-matching alleles (e.g. spike-in controls) stay in the output but get zeroed locus frequency stats."
            arg_type = String
            default = ""
        end

        add_arg_group!(ex, "RSS / core extraction", "exact_rss")
        @add_arg_table! ex begin
        "--rss"
            help = "Comma-separated list of rss fragments: heptamer, spacer, nonamer"
            default = "heptamer"
            arg_type = String
        "-a", "--affix"
            help = "Number of bases to extract from the non-RSS side of the sequence"
            arg_type = Int
            default = 13
            range_tester = (x->x >= 1)
        "--extension"
            help = "Length of extension on RSS side instead of RSS elements"
            arg_type = Int
        "--border"
            help = "Number of nucleotides from both read ends forming a border; reject if extension overlaps"
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "--adjust-per-gene-extension"
            help = "Auto-reduce extension per gene and side to avoid crossing the border"
            action = :store_true
        "--adjust-percent"
            help = "Target fraction (0-1] of reads per gene that should be safe (no border overlap) when calibrating per-gene extension"
            arg_type = Float64
            default = 1.0
            range_tester = (x-> (x > 0.0) & (x <= 1.0))
        end

        add_arg_group!(ex, "Count and allelic-ratio filters (stage 1)", "exact_filters")
        @add_arg_table! ex begin
        "-c", "--min-count"
            help = "Filter count: reads per (donor, allele, sequence), collapsing flank variants. 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-fullcount"
            help = "Filter full_count: identical rows (sequence + flanks) per donor. 0 = off. Default 5."
            default = 5
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-seqlen"
            help = "Drop candidates whose matched sequence is shorter than this (nt). 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "-f", "--min-allelic-ratio"
            help = "Filter allelic_ratio: count÷max(count in donor+gene). IgDiscover allele_ratio (÷max). 0 = off. Default 0.1."
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-e", "--expect"
            help = "Optional TSV (columns: name, ratio): per-gene/allele floor for --min-allelic-ratio."
            arg_type = String
        "--min-full-allelic-ratio"
            help = "Filter full_allelic_ratio: full_count÷max(full_count in donor+gene). 0 = off. Default 0.1."
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--expect-full"
            help = "Optional TSV (columns: name, ratio): per-gene/allele floor for --min-full-allelic-ratio."
            arg_type = String
        "--min-recurrence"
            help = "Require candidate in at least this many donors (n_donors). 0 = off. n_donors is written only with --diagnostic."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-peak-allelic-ratio"
            help = "Filter peak_allelic_ratio = max(full_allelic_ratio) across donors. 0 = off. Column written only with --diagnostic."
            default = 0.0
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-reads-total"
            help = "Filter n_reads_total: sum of full_count for this sequence across the run. 0 = off. Column written only with --diagnostic."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        end

        add_arg_group!(ex, "Gene-usage frequency filters (stage 2)", "exact_ref")
        @add_arg_table! ex begin
        "--min-gene-fraction"
            help = "÷sum filter: gene_fraction = count÷sum(accepted count in donor+gene). Not IgDiscover allelic_ratio. 0 = off. Column written only with --diagnostic."
            default = 0.0
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-d", "--deletion"
            help = "Optional TSV (columns: name, ratio): per-gene/allele floor for --min-gene-case-freq."
            arg_type = String
        "--min-gene-case-freq"
            help = "gene_case_freq = gene_count÷case_count in donor. Flags possible gene deletions. 0 = off. Column written only with --diagnostic."
            default = 0.0
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-allele-cohort-fold"
            help = "Min allele_cohort_fold (this donor's count ÷ cohort-median of this allele). 0 = off. Column written only with --diagnostic."
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "--min-gene-cohort-fold"
            help = "Min gene_cohort_fold (this donor's gene_count ÷ cohort-median). 0 = off. Column written only with --diagnostic."
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "-r", "--refgene"
            help = "Space separated reference genes to use for computing ratio"
            nargs = '*'  # Accepts zero or more values
            arg_type = String
        "--ref-fasta"
            help = "Optional reference FASTA file to check if sequences are in database (adds isin_db column)"
            arg_type = String
        end

        add_arg_group!(ex, "Output", "exact_out")
        @add_arg_table! ex begin
        "-t", "--top"
            help = "At most N flank variants per allele (1 = collapsed)."
            arg_type = Int
            default = 1
            range_tester = (x->x >= 1)
        "-l", "--limit"
            help = "Limit input reads; 0 = no limit (for testing)."
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "--raw"
            help = "Write uncollapsed per-match TSV (before count/ratio filters). Separate from --diagnostic."
            arg_type = String
        "--diagnostic"
            help = "Write intermediate statistics columns (cohort medians/folds, chimera_score, n_donors, gene_count, …) on both the filtered TSV and the .full.tsv.gz table. Default output is identifiers, counts, allelic ratios, and sequence/flanks. Filters still use the hidden columns; see reject_reason on the full table."
            action = :store_true
        "-n", "--noplot"
            help = "Disable unicode gene plot"
            action = :store_true
        end

        ex.description = "Exact match of demultiplexed reads to database alleles. Writes a slim TSV by default plus <output>.full.tsv.gz with every candidate and reject_reason."
        ex.epilog = """
Default TSV columns: well, case, gene, db_name, count, full_count, allelic_ratio, full_allelic_ratio, sequence and flanks (plus isin_db / refgene ratios when those flags are set). The full table also has reject_reason / reject_stage.

Filters (count, allelic ratio, cohort fold, …) always run; their intermediate columns are omitted unless you pass --diagnostic. --raw PATH dumps every uncollapsed match before filtering (a different, much wider table).

Examples:
  immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V
  immunediscover search exact demux.tsv.gz IGHV.fasta exact_V.tsv.gz -g V --diagnostic
"""

    return s
end
