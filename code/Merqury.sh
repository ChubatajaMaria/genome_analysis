#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00
#SBATCH -J merqury_unpolished
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Set variables
WORK_DIR="/home/mach8934/genome_analysis/data/merqury_eval"
UNPOLISHED_ASSEMBLY="/home/mach8934/genome_analysis/data/assembly.fasta"
# Source k-mer database files
SOURCE_KMER_DB="/proj/uppmax2025-3-3/Genome_Analysis/4_Zhou_2023/reads/chr3_nanopore_15mers.ktab"
# Local copy in working directory
TOOLS_SIF="/proj/uppmax2025-3-3/Genome_Analysis/tools.sif"
OUTPUT="unpolished_eval"

# Create output directory
mkdir -p $WORK_DIR
cd $WORK_DIR

# Run MerquryFK
singularity exec $TOOLS_SIF FastK -t -T2 -k15 -N"/proj/uppmax2025-3-3/GA_mach8934/merqury/nanopore_db" /home/mach8934/genome_analysis/data/raw_data/chr3_clean_nanopore.fq.gz
singularity exec $TOOLS_SIF MerquryFK "/proj/uppmax2025-3-3/GA_mach8934/merqury/nanopore_db" "/home/mach8934/genome_analysis/data/assembly.fasta" "unpolished_eval"

