#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 30:00
#SBATCH -J quast_scaffold
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=/home/mach8934/genome_analysis/data/hic_scaffolding/evaluation/%x.%j.out

# Load modules
module load bioinfo-tools
module load quast/5.0.2

# Set variables
SRC_DIR="/home/mach8934/genome_analysis/data/hic_scaffolding"
JOB_DIR="$SRC_DIR/evaluation"
SCAFFOLDS="$SRC_DIR/yahs.out_scaffolds_final.fa"

# Create the output directory if it doesn't exist
mkdir -p $JOB_DIR

# Change to the job directory
cd $JOB_DIR

# Run QUAST for assembly evaluation
echo "Running QUAST evaluation..."
python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py $SCAFFOLDS \
  --eukaryote \
  --threads 4 \
  --large \
  --labels Scaffold \
  -o $JOB_DIR

echo "QUAST evaluation completed successfully!"
