function add_discover_args!(s)
        @add_arg_table! s["discover"] begin
            "blast"
                help = "BLAST-based candidate discovery with trimming, filtering, and identity clustering"
                action = :command
            "hsmm"
                help = "Detect D genes using an HSMM fit from reference D RSS flanks (V/J masked)"
                action = :command
            "selftest"
                help = "Score recovery of known-novel alleles from a discovery full table (base vs truth FASTA)"
                action = :command
        end

        bl = s["discover"]["blast"]

        @add_arg_table! bl begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with database sequences"
            required = true
        "output"
            help = "TSV file to save discovery results"
            required = true
        end

        add_arg_group!(bl, "Inputs and outputs", "blast_io")
        @add_arg_table! bl begin
        "-p", "--pseudo"
            help = "FASTA file with pseudo-genes"
            arg_type = String
            default = ""
        "--full-output"
            help = "Path for the full annotated candidate table (every candidate + reject_reason/reject_stage). Default: <output>.full.tsv.gz"
            arg_type = String
        "--work-dir"
            help = "Directory for BLAST cache, temporary query FASTA, combined/extended DB, and affix files (relative paths use pwd()). Nothing is written beside the input TSV."
            default = blast_default("work-dir")
            arg_type = String
        end

        add_arg_group!(bl, "Gene preset", "blast_preset")
        @add_arg_table! bl begin
        "-g", "--gene"
            help = "Use gene preset parameters for V, D, or J analysis (can be overridden by explicit parameters)"
            arg_type = String
            range_tester = (x->x ∈ keys(BLAST_PRESETS))
        "-G", "--show-presets"
            help = "Show preset parameters for V, D, and J analysis"
            action = :store_true
        end

        add_arg_group!(bl, "Extension and trimming", "blast_trim")
        @add_arg_table! bl begin
        "--forward"
            help = "Forward extension length"
            default = blast_default("forward")
            arg_type = Int
            range_tester = (x->x >= 0)
        "--reverse"
            help = "Reverse extension length"
            default = blast_default("reverse")
            arg_type = Int
            range_tester = (x->x >= 0)
        "-q", "--minquality"
            help = "Minimum fraction (0–1) of affix positions that match the read in the semi-global affix–read alignment used for 5'/3' trimming; prefix and suffix each must meet this or the candidate is dropped."
            default = blast_default("minquality")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-corecov"
            help = "Minimum ratio length(aln_qseq)/length(db_seq) after trimming"
            default = blast_default("min-corecov")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        add_arg_group!(bl, "BLAST search", "blast_search")
        @add_arg_table! bl begin
        "-a", "--args"
            help = "Additional arguments to pass to blastn"
            arg_type = String
            default = blast_default("args")
        "-d", "--max-blast-mismatch"
            help = "Max BLAST mismatch per cluster (pre-trim; drops rows before the full table)"
            default = blast_default("max-blast-mismatch")
            arg_type = Int
            range_tester = (x->x >= 0)
        "--max-aln-mismatch"
            help = "Max edit distance of trimmed core (aln_qseq) vs reference allele (output filter)"
            default = blast_default("max-aln-mismatch")
            arg_type = Int
            range_tester = (x->x >= 0)
        "-e", "--edge"
            help = "Minimum number of nucleotides required between target gene and end of the read"
            default = blast_default("edge")
            arg_type = Int
            range_tester = (x->x >= 0)
        "-s", "--subjectcov"
            help = "Minimum subject (database) coverage"
            default = blast_default("subjectcov")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-read-length"
            help = "Drop reads shorter than this (nt) BEFORE BLAST (0 = off). Speeds BLAST and removes short-read noise; folded into the BLAST cache key."
            default = blast_default("min-read-length")
            arg_type = Int
            range_tester = (x -> x >= 0)
        end

        add_arg_group!(bl, "Cluster and output filters", "blast_filters")
        @add_arg_table! bl begin
        "-l", "--length"
            help = "Minimum trimmed core length (nt)."
            default = blast_default("length")
            arg_type = Int
            range_tester = (x->x >= 1)
        "-c", "--min-count"
            help = "Filter count: reads per (donor, allele, trimmed core), summing cluster variants. 0 = off."
            default = blast_default("min-count")
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-fullcount"
            help = "Filter full_count: reads in one BLAST cluster (allele + raw hit). 0 = off."
            default = blast_default("min-fullcount")
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-allelic-ratio"
            help = "Filter allelic_ratio: count÷max(count in donor+gene). 0 = off."
            default = blast_default("min-allelic-ratio")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-f", "--min-full-allelic-ratio"
            help = "Filter full_allelic_ratio: full_count÷max(full_count in donor+gene). 0 = off."
            default = blast_default("min-full-allelic-ratio")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-peak-allelic-ratio"
            help = "Filter peak_allelic_ratio = max(full_allelic_ratio) across donors. 0 = off."
            default = blast_default("min-peak-allelic-ratio")
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-reads-total"
            help = "Filter n_reads_total across the run. 0 = off."
            default = blast_default("min-reads-total")
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-recurrence"
            help = "Require candidate in at least this many donors (n_donors). 0 = off."
            default = blast_default("min-recurrence")
            arg_type = Int
            range_tester = (x -> x >= 0)
        "--max-homopolymer"
            help = "Drop candidates with homopolymer run longer than this. 0 = off."
            default = blast_default("max-homopolymer")
            arg_type = Int
            range_tester = (x -> x >= 0)
        "-i", "--isin"
            help = "On by default: a non-exact candidate whose trimmed sequence is an exact substring of a known allele is labelled with that allele. Pass -i/--isin to disable, always emitting a novel hashed name instead."
            action = :store_false
        "--keep-failed"
            help = "Keep rows where trimming failed (aln_qseq empty). By default such rows are dropped."
            action = :store_true
        end

        add_arg_group!(bl, "Run control", "blast_run")
        @add_arg_table! bl begin
        "-o", "--overwrite"
            help = "Overwrite existing files (i.e BLAST cache)"
            action = :store_true
        "-v", "--verbose"
            help = "Print verbose output and save intermediate files"
            action = :store_true
        end

        hs = s["discover"]["hsmm"]

        @add_arg_table! hs begin
        "tsv"
            help = "TSV/TSV.GZ demultiplex file with columns well, case, name, genomic_sequence"
            required = true
        "fasta"
            help = "FASTA file with D alleles (known reference)"
            required = true
        "output"
            help = "TSV.GZ file to save detected D alleles (novel and/or known) with flanks"
            required = true
        end

        add_arg_group!(hs, "Reference D selection", "hsmm_ref")
        @add_arg_table! hs begin
        "--select-min-count"
            help = "Reference D selection: minimum count per exact match (collapsed per sequence). 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--select-min-fullcount"
            help = "Reference D selection: minimum full_count per exact match row. 0 = off."
            default = 10
            arg_type = Int
            range_tester = (x->x >= 0)
        "--select-min-allelic-ratio"
            help = "Reference D selection: allelic_ratio ≥ T (count÷max in donor+gene). 0 = off."
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-l", "--limit"
            help = "Limit number of demultiplexed reads to process (0 means no limit)"
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        end

        add_arg_group!(hs, "HSMM model", "hsmm_model")
        @add_arg_table! hs begin
        "--min-posterior"
            help = "Minimum posterior probability for accepting an HSMM detection"
            default = 0.7
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-gene-len"
            help = "Minimum D gene length for HSMM duration model (auto if 0)"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 0)
        "--max-gene-len"
            help = "Maximum D gene length for HSMM duration model (auto if 0)"
            default = 70
            arg_type = Int
            range_tester = (x->x >= 0)
        end

        add_arg_group!(hs, "Output filters", "hsmm_out")
        @add_arg_table! hs begin
        "--min-count"
            help = "Output filter: count = HSMM detections clearing --min-posterior per (well, case, sequence). 0 = off."
            default = 10
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-fullcount"
            help = "Output filter: full_count = all HSMM detections per (well, case, sequence). 0 = off."
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--min-allelic-ratio"
            help = "Output filter: allelic_ratio = count÷max(count in donor+gene). 0 = off."
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-heptamer-prob-pre"
            help = "Minimum probability under pre-heptamer PWM to keep detection (0 disables)"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-heptamer-prob-post"
            help = "Minimum probability under post-heptamer PWM to keep detection (0 disables)"
            default = 0.05
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        @add_arg_table! s["discover"]["selftest"] begin
        "discovery"
            help = "Discovery FULL table (TSV/TSV.GZ) with reject_reason/reject_stage, e.g. <output>.full.tsv.gz"
            required = true
        "base"
            help = "BASE reference FASTA used for discovery (known alleles only)"
            required = true
        "truth"
            help = "TRUTH FASTA (known + novel); truth-novel = sequences not in BASE"
            required = true
        "output"
            help = "TSV path for the per-allele recovery table"
            required = true
        "--seq-col"
            help = "Discovery column holding the candidate core sequence"
            default = "aln_qseq"
            arg_type = String
        "--no-substring"
            help = "Require exact sequence match (disable substring matching)"
            action = :store_true
        "--metrics-output"
            help = "Optional TSV path to save the metric-separation table (which threshold best splits true from false novel candidates)"
            arg_type = String
        "-g", "--gene"
            help = "Gene preset used for discover blast (optional; reconstructs filter thresholds for marginal filter audit)"
            default = ""
            arg_type = String
            range_tester = (x -> isempty(x) || x ∈ keys(BLAST_PRESETS))
        end

    return s
end
