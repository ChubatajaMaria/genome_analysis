#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH -J modeling
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools RepeatModeler/2.0.4


# Your commands
BuildDatabase -name niphotrichum_japonicum /home/mach8934/genome_analysis/data/assembly.fasta \

RepeatModeler -database niphotrichum_japonicum -threads 16 -LTRStruct

