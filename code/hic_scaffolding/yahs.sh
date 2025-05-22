#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J scaffolding
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load samtools/1.20  # For creating the FASTA index

# Set variables
WORK_DIR="/proj/uppmax2025-3-3/GA_mach8934/hic_scaffolding"
REF_GENOME="/home/mach8934/genome_analysis/data/polished_data/polished_assembly.fasta"
HIC_TO_CONTIGS="/home/mach8934/genome_analysis/data/hic_scaffolding/hic_to_contigs.bam"
YAHS_PATH="/proj/uppmax2025-3-3/Genome_Analysis/yahs/yahs"

# Make sure the output directory exists
mkdir -p $WORK_DIR

# Create FASTA index file if it doesn't exist
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Creating FASTA index file..."
    samtools faidx $REF_GENOME
fi

# Check that the index was created successfully
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "ERROR: Failed to create FASTA index file. Exiting."
    exit 1
fi

# Change to the work directory
cd $WORK_DIR

# Check if the YaHS executable exists and is executable
if [ ! -x "$YAHS_PATH" ]; then
    echo "ERROR: YaHS executable not found or not executable at $YAHS_PATH"
    echo "Adding execute permission and trying again..."
    chmod +x $YAHS_PATH
fi

# Run YaHS for scaffolding
echo "Starting YaHS scaffolding..."
$YAHS_PATH -q 10 $REF_GENOME $HIC_TO_CONTIGS

# Check if YaHS completed successfully
if [ $? -ne 0 ]; then
    echo "ERROR: YaHS scaffolding failed."
    exit 1
fi

# Create symbolic links to the results
echo "Creating symbolic links to results..."
ln -sf "$WORK_DIR/yahs.out.bin" "/home/mach8934/genome_analysis/data/hic_scaffolding/yahs.out.bin"
ln -sf "$WORK_DIR/yahs.out_scaffolds_final.fa" "/home/mach8934/genome_analysis/data/hic_scaffolding/yahs.out_scaffolds_final.fa"
ln -sf "$WORK_DIR/yahs.out_scaffolds_final.agp" "/home/mach8934/genome_analysis/data/hic_scaffolding/yahs.out_scaffolds_final.agp"

echo "YaHS scaffolding completed successfully!"
