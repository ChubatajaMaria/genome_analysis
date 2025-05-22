#!/bin/bash -l
#SBATCH -A uppmax2025-3-3
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 20:00:00
#SBATCH -J polishing_with_short_read
#SBATCH --mail-type=ALL
#SBATCH --mail-user mariia.chubataia.8934@student.uu.se
#SBATCH --output=%x.%j.out
# Load modules
module load bioinfo-tools Pilon/1.24


# Your commands
java -Xmx16G -jar $PILON_HOME/pilon.jar \
--genome /home/mach8934/genome_analysis/data/assembly.fasta \
--frags /home/mach8934/genome_analysis/data/03_out_samtools/aligned.bam \
--output /home/mach8934/genome_analysis/data/03_polished_data/polished_assembly \
--threads 8 \
--changes
