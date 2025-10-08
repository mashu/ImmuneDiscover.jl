# Typical Workflows

Complete analysis pipelines for common use cases.

---

## Workflow 1: V Gene Discovery

Complete pipeline from raw reads to validated novel alleles.

```bash
# 1. Demultiplex plate
immunediscover demultiplex plate_V.fastq.gz indices.tsv demux_V.tsv.gz --length 250

# 2. Exact search for known alleles
immunediscover search exact demux_V.tsv.gz IGHV_known.fasta exact_V.tsv.gz -g V

# 3. BLAST discovery of novel alleles
immunediscover search blast demux_V.tsv.gz IGHV_known.fasta blast_V.tsv.gz -g V

# 4. Export novel alleles to FASTA
immunediscover table fasta blast_V.tsv.gz novel_V_candidates.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 2 --cleanup " Novel"

# 5. BWA QC - filter off-target sequences
bwa index GCF_000001405.25.fasta
# First create TSV for BWA input
immunediscover table select blast_V.tsv.gz bwa_input.tsv.gz \
  --columns allele_name,prefix,aln_qseq,suffix
# Map to genome
immunediscover analyze bwa bwa_input.tsv.gz bwa_validated.tsv.gz \
  GCF_000001405.25.fasta --chromosome "chromosome 14" \
  --colname allele_name --colseq prefix,aln_qseq,suffix

# 6. Export validated novel alleles
immunediscover table fasta bwa_validated.tsv.gz novel_V.fasta \
  --colname allele_name --colseq aln_qseq --cleanup " Novel"

# 7. Merge with known database
immunediscover fasta merge IGHV_updated.fasta IGHV_known.fasta novel_V.fasta

# 8. Re-run exact search with updated database
immunediscover search exact demux_V.tsv.gz IGHV_updated.fasta exact_V_final.tsv.gz -g V

# 9. Association analysis to identify haplotype blocks
immunediscover analyze association exact_V_final.tsv.gz association_edges.tsv.gz \
  --min-donors 3 --min-support 5 --clusters association_clusters.tsv

# 10. Haplotype inference
immunediscover analyze haplotype exact_V_final.tsv.gz haplotypes.tsv \
  --novel-fasta novel_V.fasta --mincount 5 --min-ratio 0.1
```

---

## Workflow 2: D Gene Discovery

D genes are short and highly variable. Two approaches available:

### Option A: BLAST with Extension

```bash
# 1. Demultiplex
immunediscover demultiplex plate_D.fastq.gz indices.tsv demux_D.tsv.gz --length 200

# 2. BLAST discovery with D gene preset
immunediscover search blast demux_D.tsv.gz IGHD_known.fasta blast_D.tsv.gz -g D

# 3. Export novel D genes
immunediscover table fasta blast_D.tsv.gz novel_D.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 3

# 4. Merge databases
immunediscover fasta merge IGHD_updated.fasta IGHD_known.fasta novel_D.fasta --cleanup " Novel"
```

### Option B: HSMM Detection (Recommended)

Better for short D genes when V/J mask D regions.

```bash
# 1. Demultiplex
immunediscover demultiplex plate_D.fastq.gz indices.tsv demux_D.tsv.gz --length 200

# 2. HSMM detection
immunediscover search hsmm demux_D.tsv.gz IGHD_known.fasta hsmm_D.tsv.gz \
  --ratio 0.2 --mincount 10 --min-posterior 0.7

# 3. Export novel D genes
immunediscover table fasta hsmm_D.tsv.gz novel_D.fasta \
  --colname allele_name --colseq sequence --filter "Novel" --mincase 3

# 4. Merge databases
immunediscover fasta merge IGHD_updated.fasta IGHD_known.fasta novel_D.fasta
```

---

## Workflow 3: Multi-Plate Batch Processing

Process multiple plates and aggregate results.

```bash
# Demultiplex all plates
for fastq in plate*.fastq.gz; do
  base=${fastq%.fastq.gz}
  immunediscover demultiplex $fastq ${base}_indices.tsv ${base}_demux.tsv.gz
done

# Exact search on all plates
for demux in *_demux.tsv.gz; do
  base=${demux%_demux.tsv.gz}
  immunediscover search exact $demux IGHV.fasta ${base}_exact.tsv.gz -g V
done

# Collect all results into single file
immunediscover table collect "plate*_exact.tsv.gz" all_plates_exact.tsv.gz

# Aggregate by allele across all plates
immunediscover table aggregate all_plates_exact.tsv.gz allele_summary.tsv.gz \
  --group-by db_name,sequence --keep-columns gene

# Filter high-confidence alleles
immunediscover table filter allele_summary.tsv.gz high_conf.tsv.gz \
  --column count --operator ">=" --threshold 100

# Export to FASTA
immunediscover table fasta high_conf.tsv.gz validated_alleles.fasta \
  --colname db_name --colseq sequence
```

---

## Workflow 4: Complete V-D-J Analysis

Analyze all gene segments from a complete library.

```bash
# 1. Demultiplex V and D/J libraries separately
immunediscover demultiplex plate_V.fastq.gz indices.tsv demux_V.tsv.gz
immunediscover demultiplex plate_DJ.fastq.gz indices.tsv demux_DJ.tsv.gz

# 2. V gene analysis
immunediscover search exact demux_V.tsv.gz IGHV.fasta exact_V.tsv.gz -g V
immunediscover search blast demux_V.tsv.gz IGHV.fasta blast_V.tsv.gz -g V

# 3. D gene analysis
immunediscover search hsmm demux_DJ.tsv.gz IGHD.fasta hsmm_D.tsv.gz

# 4. J gene analysis  
immunediscover search exact demux_DJ.tsv.gz IGHJ.fasta exact_J.tsv.gz -g J

# 5. Export novel alleles
immunediscover table fasta blast_V.tsv.gz novel_V.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 2
immunediscover table fasta hsmm_D.tsv.gz novel_D.fasta \
  --colname allele_name --colseq sequence --filter "Novel" --mincase 3

# 6. Association analysis for V genes
immunediscover analyze association exact_V.tsv.gz association_V.tsv.gz \
  --min-donors 3 --clusters clusters_V.tsv

# 7. Haplotype inference across all genes
immunediscover table collect "exact_*.tsv.gz" all_exact.tsv.gz
immunediscover analyze haplotype all_exact.tsv.gz haplotypes_all.tsv \
  --novel-fasta novel_V.fasta
```

---

## Workflow 5: Quality Control Pipeline

Comprehensive QC before finalizing novel alleles.

```bash
# 1. Initial discovery
immunediscover search blast demux.tsv.gz IGHV.fasta blast.tsv.gz -g V

# 2. Check demux statistics
cat demux.tsv.gz.log  # Review length distribution and counts per well

# 3. Filter for high-quality candidates
immunediscover table filter blast.tsv.gz filtered.tsv.gz \
  --column count --operator ">=" --threshold 10

# 4. BWA genome mapping QC
bwa index genome.fasta
immunediscover analyze bwa filtered.tsv.gz bwa_qc.tsv.gz genome.fasta

# 5. Association analysis to detect artifacts
immunediscover analyze association bwa_qc.tsv.gz association.tsv.gz \
  --min-support 5 --min-jaccard 0.3

# 6. Review association clusters for unexpected patterns
cat association_clusters.tsv

# 7. Haplotype inference to check heterozygosity rates
immunediscover analyze haplotype bwa_qc.tsv.gz haplotypes.tsv

# 8. Export validated novel alleles
immunediscover table fasta bwa_qc.tsv.gz validated_novel.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --mincase 3 --cleanup " Novel"
```

---

## Workflow 6: Comparing Databases

Compare and update allele databases.

```bash
# 1. Discover alleles from your data
immunediscover search blast demux.tsv.gz IGHV_v1.fasta blast.tsv.gz -g V
immunediscover table fasta blast.tsv.gz discovered.fasta \
  --colname allele_name --colseq aln_qseq --filter "Novel" --cleanup " Novel"

# 2. Compare with IMGT database
immunediscover fasta diff IGHV_v1.fasta IGHV_IMGT.fasta discovered.fasta

# 3. Merge unique sequences
immunediscover fasta merge IGHV_combined.fasta IGHV_v1.fasta IGHV_IMGT.fasta discovered.fasta

# 4. Add hash suffixes for version control
immunediscover fasta hash IGHV_combined.fasta > IGHV_combined_hashed.fasta

# 5. Verify with exact search
immunediscover search exact demux.tsv.gz IGHV_combined_hashed.fasta final.tsv.gz -g V
```

---

## Workflow Tips

### Iterative Discovery

For best results, iterate discovery and exact search:

```bash
# Round 1: Initial discovery
immunediscover search blast demux.tsv.gz DB_v1.fasta blast_r1.tsv.gz -g V
immunediscover table fasta blast_r1.tsv.gz novel_r1.fasta --filter "Novel" --mincase 3
immunediscover fasta merge DB_v2.fasta DB_v1.fasta novel_r1.fasta

# Round 2: Re-search with updated DB
immunediscover search exact demux.tsv.gz DB_v2.fasta exact_r2.tsv.gz -g V
immunediscover search blast demux.tsv.gz DB_v2.fasta blast_r2.tsv.gz -g V
immunediscover table fasta blast_r2.tsv.gz novel_r2.fasta --filter "Novel" --mincase 5
immunediscover fasta merge DB_v3.fasta DB_v2.fasta novel_r2.fasta

# Round 3: Final validation
immunediscover search exact demux.tsv.gz DB_v3.fasta exact_final.tsv.gz -g V
```

### Parallel Processing

Process gene types in parallel:

```bash
# Start all in background
immunediscover search exact demux_V.tsv.gz IGHV.fasta exact_V.tsv.gz -g V &
immunediscover search hsmm demux_DJ.tsv.gz IGHD.fasta hsmm_D.tsv.gz &
immunediscover search exact demux_DJ.tsv.gz IGHJ.fasta exact_J.tsv.gz -g J &

# Wait for completion
wait

# Collect results
immunediscover table collect "exact_*.tsv.gz" all_genes.tsv.gz
```

### Data Organization

Recommended file naming convention:

```
project/
  ├── raw/
  │   ├── plate1_V.fastq.gz
  │   ├── plate1_DJ.fastq.gz
  │   └── indices.tsv
  ├── demux/
  │   ├── plate1_V_demux.tsv.gz
  │   └── plate1_DJ_demux.tsv.gz
  ├── exact/
  │   ├── plate1_V_exact.tsv.gz
  │   └── plate1_DJ_exact.tsv.gz
  ├── discovery/
  │   ├── plate1_V_blast.tsv.gz
  │   └── plate1_DJ_hsmm.tsv.gz
  ├── novel/
  │   ├── novel_V.fasta
  │   └── novel_D.fasta
  ├── databases/
  │   ├── IGHV_v1.fasta
  │   ├── IGHV_v2.fasta (after iteration)
  │   └── IGHD_v1.fasta
  └── analysis/
      ├── association_edges.tsv.gz
      ├── haplotypes.tsv
      └── final_validated.fasta
```

