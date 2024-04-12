# Snakemake-ChIPseq
Snakemake workflow for ChIPseq analysis and plot replication from the published paper: Barutcu et al. RUNX1 Contributes to Higher-Order Chromatin Organization and Gene Regulation in Breast Cancer Cells. Biochimica et Biophysica Acta 1859 (11): 1389–97. PMID: PMID: 27514584

Deliverables:
  - QC and read trimming with FastQC, MultiQC, and trimmomatic
  - Build a full genome index using bowtie2-build for hg38
  - Align your reads to the human reference genome with Bowtie-2 using the generated index
  - Sort and index all generated alignment files using Samtools
  - Perform quality control with flagstats and MultiQC
  - Generate bigwig files and comparing “similarity” between samples
  - Run multiBigwig summary to generate a single matrix containing all of the information the bigwigf iles
  - Run plotCorrelation to generate a figure displaying the pearson correlation values between all of the samples
  - Perform peak calling using HOMER
  - Determine a set of reproducible peaks using bedtools
  - filter out any reproducible peaks that fall into blacklisted regions
  - Annotate peaks to their nearest genomic feature
  - Perform motif finding to look for enrichment of known binding sequences
  - Run the computeMatrix and plotProfile utilities in deeptools to identify signal coverage found in the gene bodies of all genes in the reference genome
  - Reproduce figure 2f: a stacked bar chart displaying how often a Runx1 binding peak was observed in genes that were found to be differentially expressed in the same cell line upon Runx1 knockdown using a shRNA
  - generate visualizations of the peaks found in the promoter regions of two key genes reported by the paper (MALAT1 and NEAT1)
