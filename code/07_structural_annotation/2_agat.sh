#!/bin/bash
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -J agat
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load Modules
module load bioinfo-tools AGAT/1.3.2

# Define file paths
GFF3_FILE="/home/mach8934/genome_analysis/data/07_structural_annotation/braker.gff3"
GFF_CONVERTED="/home/mach8934/genome_analysis/data/07_structural_annotation/agat/braker_converted.gff3"
GFF_STD_NAME="braker_standardized.gff3"
GFF_STANDARD_PATH="/home/mach8934/genome_analysis/data/07_structural_annotation/agat/$GFF_STD_NAME"
GENOME_FILE="/home/mach8934/genome_analysis/data/03_polished_data/polished_assembly.fasta"

# Create output directory if it doesn't exist
mkdir -p "/home/mach8934/genome_analysis/data/07_structural_annotation/agat"

# Change to output directory
cd "/home/mach8934/genome_analysis/data/07_structural_annotation/agat"

# Convert GFF into GFF (bioperl) format
agat_convert_sp_gxf2gxf.pl -gff $GFF3_FILE -o $GFF_CONVERTED

# Manage IDs
agat_sp_manage_IDs.pl -gff $GFF_CONVERTED -o $GFF_STANDARD_PATH --prefix NJAP --ensembl

# Get feature statistics
agat_sp_statistics.pl --gff $GFF_STANDARD_PATH -o braker_standardized.statistics.txt

# Extract protein sequences
agat_sp_extract_sequences.pl --gff $GFF_STANDARD_PATH --fasta $GENOME_FILE -p -o braker_standardized.aa
