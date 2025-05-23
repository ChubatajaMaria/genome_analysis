#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -J bam_coverage
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load samtools/1.20

WORK_DIR="/home/mach8934/genome_analysis/data/06_rna_mapping"
cd $WORK_DIR

echo "Processing Control_1"
samtools coverage -m Control_1_sorted.bam -o Control_1.coverage.txt
echo "Processing Control_2"
samtools coverage -m Control_2_sorted.bam -o Control_2.coverage.txt
echo "Processing Control_3"
samtools coverage -m Control_3_sorted.bam -o Control_3.coverage.txt
echo "Processing Heat_treated_42_12h_1"
samtools coverage -m Heat_treated_42_12h_1_sorted.bam -o Heat_treated_42_12h_1.coverage.txt
echo "Processing Heat_treated_42_12h_2"
samtools coverage -m Heat_treated_42_12h_2_sorted.bam -o Heat_treated_42_12h_2.coverage.txt
echo "Processing Heat_treated_42_12h_3"
samtools coverage -m Heat_treated_42_12h_3_sorted.bam -o Heat_treated_42_12h_3.coverage.txt

for file in *.coverage.txt
do
    echo -n "$file: "; grep "Mean coverage:" "$file" | \
    sed -E 's/.*Mean coverage:[[:space:]]+([0-9.]+)x/\1/' | \
    awk '{sum+=$1} END {print sum/NR}'
done
