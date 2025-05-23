#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J assembly_canu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools canu/2.2

# Set variables - use project directory with more space
WORK_DIR="/proj/uppmax2025-3-3/GA_mach8934/canu"

mkdir -p $WORK_DIR
cd $WORK_DIR

# Clean up any existing files from previous runs
rm -rf chr3_assembly.*

# Run Canu assembly
canu \
  -p chr3_assembly \
  -d $WORK_DIR \
  genomeSize=50m \
  -nanopore /home/mach8934/genome_analysis/data/raw_data/chr3_clean_nanopore.fq.gz \
  useGrid=false \
  maxThreads=16

# Create symbolic links to results
echo "Creating symbolic links to results..."
mkdir -p "/home/mach8934/genome_analysis/data/canu_assembly"

# Link the main assembly output files
ln -sf "$WORK_DIR/chr3_assembly.contigs.fasta" "/home/mach8934/genome_analysis/data/canu_assembly/chr3_assembly.contigs.fasta"
ln -sf "$WORK_DIR/chr3_assembly.unassembled.fasta" "/home/mach8934/genome_analysis/data/canu_assembly/chr3_assembly.unassembled.fasta"
ln -sf "$WORK_DIR/chr3_assembly.correctedReads.fasta.gz" "/home/mach8934/genome_analysis/data/canu_assembly/chr3_assembly.correctedReads.fasta.gz"
ln -sf "$WORK_DIR/chr3_assembly.trimmedReads.fasta.gz" "/home/mach8934/genome_analysis/data/canu_assembly/chr3_assembly.trimmedReads.fasta.gz"
ln -sf "$WORK_DIR/chr3_assembly.report" "/home/mach8934/genome_analysis/data/canu_assembly/chr3_assembly.report"
