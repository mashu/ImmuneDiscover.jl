# Troubleshooting Guide

Common issues and solutions.

---

## Demultiplexing Issues

### Low Demultiplexing Rate (<50%)

**Symptoms:** Few reads assigned to wells, many unassigned.

**Solutions:**
- ✓ Verify indices.tsv has exactly 3 columns (no extras)
- ✓ Check barcode orientation (forward=5', reverse=3')
- ✓ Inspect first few reads manually for barcode positions
- ✓ Check for adapter contamination offset

### Unbalanced Wells

**Symptoms:** Some wells have very few reads.

**Solutions:**
- Review demux log file for length distributions
- Check for failed PCR reactions or pipetting errors
- Consider excluding low-count wells from downstream analysis

---

## Discovery Issues

### Too Many Novel Alleles

**Symptoms:** Hundreds of "Novel" alleles per donor.

**Solutions:**
- ✓ Increase `--mincount` (try 10+)
- ✓ Increase `--minfullfreq` (try 0.2+)
- ✓ Check `--min-corecov` (≥0.6 recommended)
- ✓ Increase `--minquality` (try 0.8+)
- ✓ Run BWA QC to filter pseudogenes
- ✓ Review cross-case median ratios (artifacts have low ratios)

### Missing Expected Alleles

**Symptoms:** Known alleles not detected.

**Solutions:**
- ✓ Temporarily lower `--mincount` to 1 for exploration
- ✓ Check `--length` filter (may exclude short alleles)
- ✓ Verify `--locus` prefix includes target genes
- ✓ Review demux length distribution
- ✓ Check database FASTA format (GENE*ALLELE)
- ✓ Examine raw exact search output (`--raw`)

### Low BLAST Coverage

**Symptoms:** Many `corecov` values <<0.6.

**Solutions:**
- Increase extension lengths (--forward, --reverse)
- Check for poor quality reads (short or truncated)
- Verify database sequences are full-length

---

## Performance Issues

### High Memory Usage

**Solutions:**
- Process plates separately, then collect:
  ```bash
  for d in plate*_demux.tsv.gz; do
    immunediscover search exact $d DB.fasta ${d%.tsv.gz}_exact.tsv.gz -g V
  done
  immunediscover table collect "plate*_exact.tsv.gz" all.tsv.gz
  ```
- Use `--limit` for testing (e.g., 10000)
- Reduce `--top` (fewer flank variants stored)

### Slow Execution

**Solutions:**
- ✓ Set `JULIA_NUM_THREADS` to available cores
- ✓ Use `--limit` for testing before full run
- ✓ Check if BLAST is bottleneck (use `--verbose`)
- ✓ For HSMM: reduce `--max-gene-len` if possible

---

## Analysis Issues

### Unexpected Associations

**Symptoms:** Alleles from different genes show high phi coefficient.

**Solutions:**
- ✓ Check for plate effects (batch artifacts)
- ✓ Verify donors are independent (no duplicates)
- ✓ Look for contamination or barcode swapping
- ✓ Increase `--min-support` to require more co-occurrences
- ✓ Filter edges by `--min-jaccard` to reduce spurious correlations

### Low Heterozygosity

**Symptoms:** Most genotypes classified as homozygous.

**Solutions:**
- Lower `--min-ratio` threshold (try 0.05)
- Check if population is actually homozygous-heavy
- Verify sequencing depth is sufficient for minor alleles
- Review allele frequency distributions

### Too Many Uncertain Genotypes

**Symptoms:** Many genotypes with >2 alleles per gene.

**Solutions:**
- Increase `--mincount` to filter sequencing errors
- Check for contamination between samples
- Verify well demultiplexing is correct
- Consider gene duplications or pseudogenes

---

## File Format Issues

### Column Not Found Errors

**Error:** `Column 'X' not found in input file`

**Solutions:**
- ✓ Use `--case-col`, `--allele-col` to specify correct names
- ✓ Verify input is output from expected prior command
- ✓ Check TSV delimiter is tab (not space or comma)

### Empty Output Files

**Symptoms:** Output TSV has headers but no data.

**Solutions:**
- Lower filtering thresholds temporarily
- Check input file has data (`zcat file.tsv.gz | head`)
- Verify database FASTA is not empty
- Review log messages for filtering statistics

---

## Quality Control Checklist

Before finalizing results, verify:

1. ✓ **Demux stats**: Check demux log for balanced well counts
2. ✓ **Gene distribution**: Plot gene counts (should be balanced if library is good)
3. ✓ **Novel allele QC**: Run BWA mapping before validation
4. ✓ **Association patterns**: Review for unexpected correlations
5. ✓ **Heterozygosity**: Check if rates match population expectations
6. ✓ **Cross-case medians**: Identify donor-specific outliers

---

## Getting Help

**If issues persist:**
- Check [GitLab issues](https://gitlab.com/gkhlab/immunediscover.jl/-/issues)
- Review command help: `immunediscover <command> --help`
- Enable verbose output: `-v` or `--verbose` (where available)
- Save intermediate files for debugging
- Check immunediscover.log in working directory

