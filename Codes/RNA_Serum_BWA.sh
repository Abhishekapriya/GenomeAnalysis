#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J Abhi_Serum_RNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

###Mapping RNA from Serum

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands

#Index reference genome as it is large
#Map RNA from Serum to genome assembly

bwa index -p EF_rna_serum_map_ERR1797969 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797969 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_1.fastq.gz /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_2.fastq.gz | samtools sort -o output_rna_serum_map_ERR1797969.bam

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands

#Index reference genome as it is large
#Map RNA from Serum to genome assembly

bwa index -p EF_rna_serum_map_ERR1797970 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797970 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_1.fastq.gz /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_2.fastq.gz | samtools sort -o output_rna_serum_map_ERR1797970.bam

# Load modules
module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
# Your commands


bwa index -p EF_rna_serum_map_ERR1797971 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797971 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_1.fastq.gz  /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_2.fastq.gz | samtools sort -o output_rna_serum_map_ERR1797971.bam
