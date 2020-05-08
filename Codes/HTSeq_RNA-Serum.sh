#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J Abhi_Serum_RNA_HTSeq
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

###Counts of Mapped Serum RNA from Serum

# Load modules

module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797969.bam | htseq-count --format bam --type CDS  /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797969.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_Serum_69_counts.txt

module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797970.bam | htseq-count --format bam --type CDS /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797970.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_Serum_70_counts.txt


module load bioinfo-tools
module load samtools/1.9
module load htseq/0.9.1

samtools sort /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797971.bam | htseq-count --format bam --type CDS /home/abhis/TranscriptomeAnalysis/HtSeq/output_rna_serum_map_ERR1797971.bam /home/abhis/GenomeAnalysis/Prokka/PROKKA_04182020/output.gff > RNA_Serum_71_counts.txt
