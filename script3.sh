#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J Abhi_Trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

# Load modules
module load bioinfo-tools
module load trimmomatic

# Your commands
trimmomatic PE -phred33 /home/abhis/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz /home/abhis/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz E745-1.L500_SZAXPI015146-56_1_clean_trim_f_paired.fq.gz E745-1.L500_SZAXPI015146-56_1_clean_trim_f_unpaired.fq.gz E745-1.L500_SZAXPI015146-56_2_clean_trim_r_paired.fq.gz E745-1.L500_SZAXPI015146-56_2_clean_trim_r_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:27 MINLEN:36
