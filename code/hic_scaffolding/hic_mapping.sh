#!/bin/bash -l

#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J hic_mapping
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out
#SBATCH --output=/home/mach8934/genome_analysis/code/hic_scaffolding/%x.%j.out

module load bioinfo-tools
module load samtools/1.20
module load bwa/0.7.18

THREADs=4

SRC_DIR="/home/mach8934/genome_analysis"
WORK_DIR="/proj/uppmax2025-3-3/GA_mach8934/hic_scaffolding"

REF_GENOME="/home/mach8934/genome_analysis/data/polished_data/polished_assembly.fasta"
FORWARD_READ="/home/mach8934/genome_analysis/data/raw_data/chr3_hiC_R1.fastq.gz"
REVERSE_READ="/home/mach8934/genome_analysis/data/raw_data/chr3_hiC_R2.fastq.gz"

cd $WORK_DIR

bwa index $REF_GENOME
bwa mem -t $THREADs $REF_GENOME $FORWARD_READ $REVERSE_READ | \
samtools view -Sb - > aln.bam

samtools collate -@ $THREADs -O -u aln.bam | \
samtools fixmate -@ $THREADs -m -u - - | \
samtools sort -@ $THREADs -u - | \
samtools markdup -@ $THREADs - $WORK_DIR/hic_to_contigs.bam
samtools index $WORK_DIR/hic_to_contigs.bam -o hic_to_contigs.bam.bai

ln -sf "$WORK_DIR/hic_to_contigs.bam" "/home/mach8934/genome_analysis/data/hic_scaffolding/hic_to_contigs.bam"
ln -sf "$WORK_DIR/hic_to_contigs.bam.bai" "/home/mach8934/genome_analysis/data/hic_scaffolding/hic_to_contigs.bam.bai"
