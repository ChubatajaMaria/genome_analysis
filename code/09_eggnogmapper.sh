#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -J eggNOGMapper
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load Modules
module load bioinfo-tools eggNOG-mapper/2.1.9

PROTEINS_FILE="/home/mach8934/genome_analysis/data/07_structural_annotation/agat/braker_standardized.aa"
GENOME_FILE="/home/mach8934/genome_analysis/data/03_polished_data/polished_assembly.fasta"

emapper.py -i "$PROTEINS_FILE" \
  --itype proteins \
  -m diamond \
  --cpu 8 \
  --go_evidence experimental \
  --output n_japonicum \
  --output_dir /home/mach8934/genome_analysis/data/08_functional_annotation\
  --decorate_gff /home/mach8934/genome_analysis/data/07_structural_annotation/agat/braker_standardized.gff3 \
  --decorate_gff_ID_field ID \
  --excel
