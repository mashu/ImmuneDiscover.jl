function add_search_args!(s)
        @add_arg_table! s["search"] begin
            "exact"
                help = "Exact match search of reads to database alleles with robust filters"
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
            help = "TSV file to save ouput"
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
            help = "Locus to filter genes to start with this string (e.g. IGHV) excluding other genes from the analysis (i.e control genes)"
            arg_type = String
            default = "IG"
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

        add_arg_group!(ex, "Count and ratio filters", "exact_filters")
        @add_arg_table! ex begin
        "-c", "--mincount"
            help = "Minimum cluster size"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
        "-f", "--minratio"
            help = "Minimum allelic ratio applied within each gene group"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-allele-mratio"
            help = "Minimum allelic ratio applied within each gene group for the allele against median"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "--min-gene-mratio"
            help = "Minimum allelic ratio applied within each gene group for the gene against median"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0))
        "--min-recurrence"
            help = "Quality filter: require a candidate sequence to appear in at least this many donors (n_donors). 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-seqlen"
            help = "Quality filter: drop candidates whose matched sequence is shorter than this (nt). 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-peak-ratio"
            help = "Quality filter: require the peak per-donor allelic ratio (max over donors of count/max-in-gene). 0 = off. Complements --minratio (per-group) with a cross-donor peak floor."
            default = 0.0
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        add_arg_group!(ex, "Reference frequency thresholds", "exact_ref")
        @add_arg_table! ex begin
        "-e", "--expect"
            help = "TSV file containing gene names and their corresponding allele_freq threshold, with two columns: name and ratio"
            arg_type = String
        "-d", "--deletion"
            help = "TSV file containing gene names and their corresponding gene_case_freq threshold, with two columns: name and ratio"
            arg_type = String
        "-r", "--refgene"
            help = "Space separated reference genes to use for computing ratio"
            nargs = '*'  # Accepts zero or more values
            arg_type = String
        "--ref-fasta"
            help = "Optional reference FASTA file to check if sequences are in database (adds isin_db column)"
            arg_type = String
        end

        add_arg_group!(ex, "Output and diagnostics", "exact_out")
        @add_arg_table! ex begin
        "-t", "--top"
            help = "Saves at most N records of flank and sequence."
            arg_type = Int
            default = 1
            range_tester = (x->x >= 1)
        "-l", "--limit"
            help = "Limit to this number of sequences, zero means no limit"
            arg_type = Int
            default = 0
            range_tester = (x->x >= 0)
        "--raw"
            help = "Unfiltered exact search results for diagnostics"
            arg_type = String
        "-n", "--noplot"
            help = "Disable unicode gene plot"
            action = :store_true
        end

    return s
end
