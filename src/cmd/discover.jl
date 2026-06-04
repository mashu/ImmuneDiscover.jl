function add_discover_args!(s)
        @add_arg_table! s["discover"] begin
            "blast"
                help = "BLAST-based candidate discovery with trimming, filtering, and identity clustering"
                action = :command
            "hsmm"
                help = "Detect D genes using an HSMM trained on RSS flanks (V/J masked)"
                action = :command
        end

        @add_arg_table! s["discover"]["blast"] begin
        "input"
            help = "TSV file with demultiplex data"
            required = true
        "fasta"
            help = "FASTA file with database sequences"
            required = true
        "output"
            help = "TSV file to save discovery results"
            required = true
        "-p", "--pseudo"
            help = "FASTA file with pseudo-genes"
            arg_type = String
            default = ""
        "-c", "--minfullcount"
            help = "Minimum full cluster size"
            default = 5
            arg_type = Int
        "-f", "--minfullratio"
            help = "Minimum allelic ratio within each gene group (count / max in gene)"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-s", "--subjectcov"
            help = "Minimum subject (database) coverage"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-d", "--maxdist"
            help = "Maximum distance allowed for alleles"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "-l", "--length"
            help = "Minimum length of the trimmed read"
            default = 290
            arg_type = Int
            range_tester = (x->x >= 1)
        "-e", "--edge"
            help = "Minimum number of nucleotides required between target gene and end of the read"
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "-a", "--args"
            help = "Additional arguments to pass to blastn"
            arg_type = String
            default = "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
        "-o",  "--overwrite"
            help = "Overwrite existing files (i.e BLAST cache)"
            action = :store_true
        "-g", "--gene"
            help = "Use gene preset parameters for V, D, or J analysis (can be overridden by explicit parameters)"
            arg_type = String
            range_tester = (x->x ∈ keys(BLAST_PRESETS))
        "-G", "--show-presets"
            help = "Show preset parameters for V, D, and J analysis"
            action = :store_true
        "--forward"
            help = "Forward extension length"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "--reverse"
            help = "Reverse extension length"
            default = 20
            arg_type = Int
            range_tester = (x->x >= 0)
        "-q", "--minquality"
            help = "Minimum fraction (0–1) of affix positions that match the read in the semi-global affix–read alignment used for 5'/3' trimming; prefix and suffix each must meet this or the candidate is dropped."
            default = 0.75
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-i", "--isin"
            help = "On by default: a non-exact candidate whose trimmed sequence is an exact substring of a known allele is labelled with that allele. Pass -i/--isin to disable, always emitting a novel hashed name instead."
            action = :store_false
        "--keep-failed"
            help = "Keep rows where trimming failed (aln_qseq empty). By default such rows are dropped."
            action = :store_true
        "--min-corecov"
            help = "Minimum ratio length(aln_qseq)/length(db_seq) after trimming"
            default = 0.6
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-v", "--verbose"
            help = "Print verbose output and save intermediate files"
            action = :store_true
        "--work-dir"
            help = "Directory for BLAST cache, temporary query FASTA, combined/extended DB, and affix files (relative paths use pwd()). Nothing is written beside the input TSV."
            default = ".immunediscover"
            arg_type = String
        end

        @add_arg_table! s["discover"]["hsmm"] begin
        "tsv"
            help = "TSV/TSV.GZ demultiplex file with columns well, case, name, genomic_sequence"
            required = true
        "fasta"
            help = "FASTA file with D alleles (known reference)"
            required = true
        "output"
            help = "TSV.GZ file to save detected D alleles (novel and/or known) with flanks"
            required = true
        "-r", "--ratio"
            help = "Allelic ratio threshold for known D selection per donor and gene"
            default = 0.2
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "-c", "--mincount"
            help = "Minimum count for a known D allele to be considered in training"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
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
        "-l", "--limit"
            help = "Limit number of demultiplexed reads to process (0 means no limit)"
            default = 0
            arg_type = Int
            range_tester = (x->x >= 0)
        "--out-mincount"
            help = "Minimum count for an extracted D (after HSMM) to keep in output"
            default = 10
            arg_type = Int
            range_tester = (x->x >= 1)
        "--out-minratio"
            help = "Minimum allelic ratio within gene (per donor) for output filtering"
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

        # Table → fasta export
    return s
end
