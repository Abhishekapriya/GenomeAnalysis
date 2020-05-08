#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J Abhi_BH_RNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

###Mapping RNA from BH

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands

#Index reference genome as it is large
#Map RNA from Serum to genome assembly

bwa index -p EF_rna_bh_map_ERR1797972 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_bh_map_ERR1797972 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_1.fastq.gz  /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_2.fastq.gz | samtools sort -O BAM -o rna_bh_map72.bam

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands

#Index reference genome as it is large
#Map RNA from Serum to genome assembly

bwa index -p EF_rna_bh_map_ERR1797973 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_bh_map_ERR1797973 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_1.fastq.gz  /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_2.fastq.gz | samtools sort -O BAM -o rna_bh_map73.bam

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands


bwa index -p EF_rna_bh_map_ERR1797974 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_bh_map_ERR1797974 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_1.fastq.gz  /home/abhis/raw_data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_2.fastq.gz | samtools sort -O BAM -o rna_bh_map74.bam
