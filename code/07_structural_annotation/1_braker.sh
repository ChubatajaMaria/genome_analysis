#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 03:00:00
#SBATCH -J braker2
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools braker/2.1.6

# Set environment variables
export AUGUSTUS_CONFIG_PATH=/home/mach8934/augustus_config
export AUGUSTUS_BIN_PATH=/sw/bioinfo/augustus/3.4.0/snowy/bin
export AUGUSTUS_SCRIPTS_PATH=/sw/bioinfo/augustus/3.4.0/snowy/scripts
export GENEMARK_PATH=/sw/bioinfo/GeneMark/4.68-es/snowy

BAM_PATH="/home/mach8934/genome_analysis/data/06_rna_mapping"
BAM_FILES=$(echo ${BAM_PATH}/*.bam | sed 's/ /,/g')

# Copy GeneMark key to home directory
cp -vf /sw/bioinfo/GeneMark/4.68-es/snowy/gm_key $HOME/.gm_key
chmod 600 $HOME/.gm_key

# Create a new working directory with a unique name
WORKING_DIR="/home/mach8934/genome_analysis/data/07_structural_annotation"
mkdir -p $WORKING_DIR

# Your commands
braker.pl \
--genome="/home/mach8934/genome_analysis/data/05_masking/out_masking/polished_assembly.fasta.masked" \
--bam="$BAM_FILES" \
--prot_seq="/home/mach8934/genome_analysis/data/proteins/embryophyte_proteomes.faa" \
--softmasking \
--etpmode \
--species="niphotrichum_japonicum" \
--useexisting \
--workingdir="$WORKING_DIR" \
--gff3 \
--cores=8
