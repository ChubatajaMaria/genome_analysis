#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J assembly_20250405
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out


#Loading of modules

module load bioinfo-tools Flye/2.9.5

#Running of gene assembly

flye --nano-raw /home/mach8934/genome_analysis/data/raw_data/chr3_clean_nanopore.fq.gz --out-dir /home/mach8934/genome_analysis/data/ --threads 8
