l#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 03:00:00
#SBATCH -J masking
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools RepeatMasker/4.1.5


# Your commands
RepeatMasker \
	-lib /home/mach8934/genome_analysis/data/05_masking/out_modeling/niphotrichum_japonicum-families.fa\
	-pa 8\
	-xsmall\
	-no_is\
	-html\
	-gff\
	-dir /home/mach8934/genome_analysis/data/05_masking/out_masking \
	/home/mach8934/genome_analysis/data/03_polished_data/polished_assembly.fasta
