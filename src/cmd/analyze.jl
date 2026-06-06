function add_analyze_args!(s)
        @add_arg_table! s["analyze"] begin
            "cooccurrence"
                help = "Analyze cross-donor co-occurrence among alleles across cases"
                action = :command
            "haplotype"
                help = "Infer approximate haplotypes per case using diploid assumptions"
                action = :command
        end

        co = s["analyze"]["cooccurrence"]
        @add_arg_table! co begin
        "input"
            help = "TSV/TSV.GZ with columns: case and db_name (or specify columns)"
            required = true
        end

        add_arg_group!(co, "Column names", "cooc_cols")
        @add_arg_table! co begin
        "-C", "--case-col"
            help = "Name of column with donor/case id"
            default = "case"
            arg_type = String
        "-A", "--allele-col"
            help = "Name of column with allele name"
            default = "db_name"
            arg_type = String
        end

        add_arg_group!(co, "Allele inclusion", "cooc_incl")
        @add_arg_table! co begin
        "-m", "--min-donors"
            help = "Minimum donors required to include an allele"
            default = 2
            arg_type = Int
            range_tester = (x->x >= 1)
        end

        add_arg_group!(co, "Clustering", "cooc_cluster")
        @add_arg_table! co begin
        "--cluster-method"
            help = "Clustering method on rho: components, complete, average, or single"
            default = "components"
            arg_type = String
            range_tester = (x-> (x ∈ ["components","complete","average","single"]))
        "--cluster-threshold"
            help = "Similarity threshold for complete-linkage clustering on rho (0..1)"
            default = 0.5
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        "--min-cluster-size"
            help = "Minimum cluster size to output"
            default = 3
            arg_type = Int
            range_tester = (x->x >= 1)
        "--clusters"
            help = "Optional path to save clusters (TSV)"
            arg_type = String
        "--debug-triangles"
            help = "Print rho-based triangle diagnostics at cluster-threshold"
            action = :store_true
        end

        hp = s["analyze"]["haplotype"]
        @add_arg_table! hp begin
        "input"
            help = "TSV/TSV.GZ with columns: case and allele (or specify columns)"
            required = true
        "output"
            help = "TSV file to save haplotype inference results"
            required = true
        end

        add_arg_group!(hp, "Column names", "hap_cols")
        @add_arg_table! hp begin
        "-C", "--case-col"
            help = "Name of column with donor/case id"
            default = "case"
            arg_type = String
        "-A", "--allele-col"
            help = "Name of column with allele name"
            default = "db_name"
            arg_type = String
        "-G", "--gene-col"
            help = "Name of column with gene name (for grouping alleles by gene)"
            default = "gene"
            arg_type = String
        end

        add_arg_group!(hp, "Genotype thresholds", "hap_thresh")
        @add_arg_table! hp begin
        "-c", "--mincount"
            help = "Minimum count for an allele to be considered"
            default = 5
            arg_type = Int
            range_tester = (x->x >= 1)
        "-r", "--min-ratio"
            help = "Minimum allelic ratio threshold (minor/major allele ratio)"
            default = 0.1
            arg_type = Float64
            range_tester = (x-> (x >= 0.0) & (x <= 1.0))
        end

        add_arg_group!(hp, "Output", "hap_out")
        @add_arg_table! hp begin
        "-f", "--novel-fasta"
            help = "Optional FASTA file with novel alleles to mark in results"
            arg_type = String
        end

    return s
end
