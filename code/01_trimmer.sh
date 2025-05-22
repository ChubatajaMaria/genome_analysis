#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -J illumina_trimmomatic
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se 
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools trimmomatic/0.39
# Your commands
trimmomatic SE -threads 4 /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R1.fastq.gz /home/mach8934/genome_analysis/data/01_trimmed_data/chr3_illumina_R1_trimmed.fq.gz HEADCROP:15 MINLEN:136
trimmomatic SE -threads 4 /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R2.fastq.gz /home/mach8934/genome_analysis/data/01_trimmed_data/chr3_illumina_R2_trimmed.fq.gz HEADCROP:15 MINLEN:136
