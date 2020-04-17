#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -J Abhi_Spades
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

# Load modules
module load bioinfo-tools
module load spades

# Your commands
spades.py -t 4 -1 /home/abhis/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz -2 /home/abhis/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz --nanopore /home/abhis/raw_data/Nanopore/E745_all.fasta.gz --careful --cov-cutoff auto -o /home/abhis/GenomeAnalysis/Spades/spades_outdir 