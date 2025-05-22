#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 03:00:00
#SBATCH -J mapping
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools RepeatMasker/
# Your commands

bwa index /home/mach8934/genome_analysis/data/assembly.fasta

bwa mem -t 8 /home/mach8934/genome_analysis/data/assembly.fasta /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R1.fastq.gz /home/mach8934/genome_analysis/data/raw_data/chr3_illumina_R2.fastq.gz | samtools sort -@ 8 -o /home/mach8934/genome_analysis/data/out_samtools/aligned.bam
samtools index /home/mach8934/genome_analysis/data/out_samtools/aligned.bam
