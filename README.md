# 🧬 RNA-Seq preprocessing Processing bash Pipeline 
This bash script automates a complete RNA-seq data preprocessing and alignment workflow for E. coli datasets using high-performance computing (HPC) with PBS job scheduling.
## Quality Control
Launches MultiQC to summarize initial data quality.
## Adapter Trimming
Uses fastp to trim sequencing adapters and improve read quality.
## Index Building
Builds a HISAT2 genome index from the E. coli reference genome.
## Read Alignment
Aligns trimmed reads to the reference genome using HISAT2 and converts .sam files to .bam with samtools.
## File Organization
Moves all resulting .bam files to a dedicated output directory.
## Gene Expression Quantification
Runs featureCounts to generate a gene count matrix from the aligned .bam files.
## Final Quality Summary
Uses MultiQC to summarize read assignment statistics from the count summary file.
