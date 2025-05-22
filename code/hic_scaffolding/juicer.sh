#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J juicer
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=/home/mach8934/genome_analysis/data/hic_scaffolding/juicer/%x.%j.out

# Load modules
module load bioinfo-tools
module load samtools/1.20

# Set variables
THREADS=8
SRC_DIR="/home/mach8934/genome_analysis/data/hic_scaffolding"
WORK_DIR="/proj/uppmax2025-3-3/GA_mach8934/hic_scaffolding/juicer"
REF_GENOME="/home/mach8934/genome_analysis/data/polished_data/polished_assembly.fasta"
BIN_FILE="$SRC_DIR/yahs.out.bin"
SCAFFOLDS_FINAL="$SRC_DIR/yahs.out_scaffolds_final.agp"

# Make sure the output directory exists
mkdir -p $WORK_DIR
mkdir -p $SRC_DIR/juicer

# Change to the work directory
cd $WORK_DIR

# Copy reference genome to work directory if needed
cp $REF_GENOME $WORK_DIR/polished_assembly.fasta
REF_GENOME="$WORK_DIR/polished_assembly.fasta"

# Step 1: Index the FASTA (needed for AGP parsing and size extraction)
echo "Indexing reference genome..."
samtools faidx $REF_GENOME

# Step 2: Create sorted alignment file for juicer_tools input
echo "Running Juicer pre command..."
/proj/uppmax2025-3-3/Genome_Analysis/yahs/juicer pre \
    $BIN_FILE $SCAFFOLDS_FINAL polished_assembly.fasta.fai | \
    sort -k2,2d -k6,6d -T ./ --parallel=$THREADS -S32G | \
    awk 'NF' > alignments_sorted.txt

# Step 3: Extract scaffold sizes from the FASTA index
echo "Extracting scaffold sizes..."
samtools faidx "$SRC_DIR/yahs.out_scaffolds_final.fa" -o "$WORK_DIR/yahs.out_scaffolds_final.fa.fai"
cut -f1,2 "$WORK_DIR/yahs.out_scaffolds_final.fa.fai" > scaffolds_final.chrom.sizes

# Step 4: Generate the .hic file for Juicebox visualization
echo "Generating .hic file..."
java -Xmx32G -jar /proj/uppmax2025-3-3/Genome_Analysis/yahs/juicer_tools_1.22.01.jar pre \
    alignments_sorted.txt yahs_output.hic scaffolds_final.chrom.sizes

# Create symbolic links to results
echo "Creating symbolic links to results..."
ln -sf "$WORK_DIR/yahs_output.hic" "$SRC_DIR/juicer/yahs_output.hic"
ln -sf "$WORK_DIR/scaffolds_final.chrom.sizes" "$SRC_DIR/juicer/scaffolds_final.chrom.sizes"
ln -sf "$WORK_DIR/alignments_sorted.txt" "$SRC_DIR/juicer/alignments_sorted.txt"

echo "Juicer processing completed successfully!"
