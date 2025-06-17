#!/usr/bin/bash
#PBS -N PL
#PBS -l nodes=1:ppn=14
#PBS -l walltime=30:00:00
#PBS -q normal
#PBS -o /nlustre/users/caivil/MASTERS/STD/qc_stdout
#PBS -e /nlustre/users/caivil/MASTERS/STD/qc_stderr
#PBS -k oe
#PBS -m ae
#PBS -M Biocaivil@gmail.com

# First run multiqc to see the quality of the data

#####################################################################
# Running Fastp for trimming adapters 
module load fastp-0.23.2

wkd="/nlustre/users/caivil/MASTERS/data/fastq_files"
cd ${wkd}

samples_1=(
  "SRR6671736" "SRR6671737" "SRR6671738" "SRR6671739" "SRR6671740"
  "SRR6671741" "SRR6671742" "SRR6671743" "SRR6671744" "SRR6671745"
  "SRR6671746" "SRR6671747" "SRR6671748" "SRR6671749" "SRR6671750"
  "SRR6671751" "SRR6671752" "SRR6671753" "SRR6671754" "SRR6671755"
  "SRR6671756" "SRR6671757" "SRR6671758" "SRR6671759" "SRR6671760"
  "SRR6671761" "SRR6671762" "SRR6671763" "SRR6671764" "SRR6671765"
  "SRR6671766" "SRR6671767" "SRR6671768" "SRR6671769" "SRR6671770"
  "SRR6671771" "SRR6671772" "SRR6671773" "SRR6671774" "SRR6671775"
  "SRR6671776" "SRR6671777" "SRR6671778" "SRR6671779" "SRR6671780"
  "SRR6671781" "SRR6671782" "SRR6671783" "SRR6671784" "SRR6671785"
  "SRR6671786" "SRR6671787" "SRR6671788" "SRR6671789" "SRR6671790"
  "SRR6671791" "SRR6671792"
)

for sample in "${samples_1[@]}"; do
  fastp -i "${sample}_1.fastq" -I "${sample}_2.fastq" \
        -o "${sample}_1_trimmed.fastq" -O "${sample}_2_trimmed.fastq"
done

##################################################################################
#Building a Hisat2 index

module load hisat2-2.1.0

wkd="/nlustre/users/caivil/MASTERS/REF/orig_index"
cd ${wkd}

hisat2-build -p 14 -f  /nlustre/users/caivil/MASTERS/REF/GCF_000750555.1_ASM75055v1_genomic.fna Ecoli_reference_hisat2

#################################################################################
#Running Hisat2 and converting .sam to .bam

module load hisat2-2.1.0
module load samtools-1.18

wkd="/nlustre/users/caivil/MASTERS/data/fastq_files"
cd ${wkd}

samples=(
  "SRR6671736"
  "SRR6671737"
  "SRR6671738"
  "SRR6671739"
  "SRR6671740"
  "SRR6671741"
  "SRR6671742"
  "SRR6671743"
  "SRR6671744"
  "SRR6671745"
  "SRR6671746"
  "SRR6671747"
  "SRR6671748"
  "SRR6671749"
  "SRR6671750"
  "SRR6671751"
  "SRR6671752"
  "SRR6671753"
  "SRR6671754"
  "SRR6671755"
  "SRR6671756"
  "SRR6671757"
  "SRR6671758"
  "SRR6671759"
  "SRR6671760"
  "SRR6671761"
  "SRR6671762"
  "SRR6671763"
  "SRR6671764"
  "SRR6671765"
  "SRR6671766"
  "SRR6671767"
  "SRR6671768"
  "SRR6671769"
  "SRR6671770"
  "SRR6671771"
  "SRR6671772"
  "SRR6671773"
  "SRR6671774"
  "SRR6671775"
  "SRR6671776"
  "SRR6671777"
  "SRR6671778"
  "SRR6671779"
  "SRR6671780"
  "SRR6671781"
  "SRR6671782"
  "SRR6671783"
  "SRR6671784"
  "SRR6671785"
  "SRR6671786"
  "SRR6671787"
  "SRR6671788"
  "SRR6671789"
  "SRR6671790"
  "SRR6671791"
  "SRR6671792"
)


for sample in "${samples[@]}"; do
  hisat2 -p 14 -x /nlustre/users/caivil/MASTERS/REF/orig_index/Ecoli_reference_hisat2 \
         -1 "${sample}_1_trimmed.fastq" -2 "${sample}_2_trimmed.fastq" -S "${sample}.sam"
  samtools view -b "${sample}.sam" > "${sample}.bam"
done
#############################################################################
# moving the output .bam files to the proper file

mv *.bam /nlustre/users/caivil/MASTERS/data/bam_files3

#########################################################################
#Running Featurecounts for counts_matrix

module load subread-2.0.8

wkd="/nlustre/users/caivil/MASTERS/data/ecoli_gft_file"
cd ${wkd}

featureCounts -p --countReadPairs -a genomic_output.saf -F SAF -o counts_x.txt -T 14 /nlustre/users/caivil/MASTERS/data/bam_files3/*.bam

############################################################################
#Running Multiqc to view assigned reads

module load multiqc
cd /nlustre/users/caivil/MASTERS/data/ecoli_gft_file
multiqc counts_x.txt.summary
