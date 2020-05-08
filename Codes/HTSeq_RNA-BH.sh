#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J Abhi_BH_RNA_HTSeq
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

###Counts of Mapped Serum RNA from Serum

# Load modules

module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map72.bam | htseq-count --format bam --type CDS  /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map72.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_BH_72_counts.txt

module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map73.bam | htseq-count --format bam --type CDS /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map73.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_BH_73_counts.txt


module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map74.bam | htseq-count --format bam --type CDS /home/abhis/TranscriptomeAnalysis/HtSeq/rna_bh_map74.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_BH_74_counts.txt