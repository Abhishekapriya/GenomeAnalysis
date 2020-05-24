#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 08:00:00
#SBATCH -J Abhi_Serum_RNA
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

###Mapping RNA from Serum

# Load modules
# Your commands

#Index reference genome as it is large
#Map RNA from Serum to genome assembly
#Count RNA from mapped serum

module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load htseq/0.9.1

bwa index -p EF_rna_serum_map_ERR1797969 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797969 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_1.fastq.gz /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_2.fastq.gz | samtools sort -O BAM -o rna_serum_map69.bam
samtools view rna_serum_map69.bam | htseq-count -f bam -t CDS rna_serum_map69.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gtf > RNA_Serum_69_counts.txt


module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load htseq/0.9.1

bwa index -p EF_rna_serum_map_ERR1797970 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797970 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_1.fastq.gz /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_2.fastq.gz | samtools sort -O BAM -o rna_serum_map70.bam
samtools view rna_serum_map70.bam | htseq-count -f bam -t CDS rna_serum_map70.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gtf > RNA_Serum_70_counts.txt


module load bioinfo-tools
module load bwa/0.7.17
module load samtools/1.9
module load htseq/0.9.1

bwa index -p EF_rna_serum_map_ERR1797971 /home/abhis/GenomeAnalysis/BWA/CombinedAlign.fasta
bwa mem -M EF_rna_serum_map_ERR1797971 /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_1.fastq.gz  /home/abhis/raw_data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_2.fastq.gz | samtools sort -O BAM -o rna_serum_map71.bam
samtools view rna_serum_map71.bam | htseq-count -f bam -t CDS rna_serum_map71.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gtf > RNA_Serum_71_counts.txt