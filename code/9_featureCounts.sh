#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH -J read_counting
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools subread/2.0.3

WORK_DIR="/home/mach8934/genome_analysis/data/09_read_counting"

# Input Files
GFF_FILE="/home/mach8934/genome_analysis/data/08_functional_annotation/n_japonicum.emapper.decorated.gff"
BAM_PATH="/home/mach8934/genome_analysis/data/06_rna_mapping"
BAM_FILES=${BAM_PATH}/*.bam

cd "$WORK_DIR"

featureCounts -p --countReadPairs \
  -T 4 \
  -a "$GFF_FILE" -t gene -g ID \
  -o "${WORK_DIR}/counts.txt" \
  $BAM_FILES
