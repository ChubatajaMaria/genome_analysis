#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 00:20:00
#SBATCH -J fastqc_20250707
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

#Load modules
module load bioinfo-tools FastQC/0.11.9

fastqc -t 4 -o /home/mach8934/genome_analysis/data/01_fastqc_data /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R1.fastq.gz
fastqc -t 4 -o /home/mach8934/genome_analysis/data/01_fastqc_data /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R2.fastq.gz
