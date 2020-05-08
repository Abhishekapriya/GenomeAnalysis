#!/bin/bash -l
#SBATCH -A g2020008
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J Abhi_Canu
#SBATCH --mail-type=ALL
#SBATCH --mail-user gau.abhi213@gmail.com

# Load modules
module load bioinfo-tools
module load canu

# Your commands
canu -p EF -d EF_Canu stopOnReadQuality=false genomeSize=3.2m -pacbio-raw /home/abhis/raw_data/PacBio/*.subreads.fastq.gz