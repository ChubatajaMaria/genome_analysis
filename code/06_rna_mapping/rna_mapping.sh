#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 09:00:00
#SBATCH -J rna_mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=/home/mach8934/genome_analysis/data/rna_mapping/%x.%j.out

# Load modules
module load bioinfo-tools samtools/1.20 HISAT2/2.2.1

# Your commands
mkdir -p /home/mach8934/genome_analysis/data/raw_data
cd /home/mach8934/genome_analysis/data/raw_data

# Build the Index
hisat2-build /home/mach8934/genome_analysis/data/03_polished_data/polished_assembly.fasta assembly_index

for FWD in /home/mach8934/genome_analysis/data/raw_data/*_f1.fq.gz; do
    SAMPLE=$(basename "$FWD" _f1.fq.gz)
    REV="/home/mach8934/genome_analysis/data/raw_data/${SAMPLE}_r2.fq.gz"
    
    hisat2 -x assembly_index -1 "$FWD" -2 "$REV" --threads 8 \
    | samtools view -bS - \
    | samtools sort -@ 8 -o "/home/mach8934/genome_analysis/data/06_rna_mapping/${SAMPLE}_sorted.bam"
    
    samtools index "/home/mach8934/genome_analysis/data/06_rna_mapping/${SAMPLE}_sorted.bam"
done
