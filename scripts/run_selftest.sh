scripts/run.sh discover blast -g V data/KI_IMD_IGHV_old38-42.demultiplexed.combined.tsv.gz data/KI+1KGP-IGHV-Base.fasta KI_IMD_IGHV.tsv.gz
scripts/run.sh discover selftest KI_IMD_IGHV.full.tsv.gz data/KI+1KGP-IGHV-Base.fasta data/KI+1KGP-IGHV-Complete.fasta --metrics-output IGHV-metrics.tsv IGHV-seltest.tsv
scripts/run.sh discover blast -g D data/KI_IMD_IGHDJ_old38-42.demultiplexed.combined.tsv.gz data/KI+1KGP-IGHD-Base.fasta KI_IMD_IGHD.tsv.gz
scripts/run.sh discover selftest KI_IMD_IGHD.full.tsv.gz data/KI+1KGP-IGHD-Base.fasta data/KI+1KGP-IGHD-Complete.fasta --metrics-output IGHD-metrics.tsv IGHD-seltest.tsv

