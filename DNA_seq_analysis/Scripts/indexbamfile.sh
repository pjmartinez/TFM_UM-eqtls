#!/bin/bash 
#SBATCH --job-name=deduplication
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH --qos=himem
#SBATCH --partition=himem
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load samtools

# "*mkdup.bam" will refer to each of the 
for file in ../align_stepwise/*_mkdup.bam
	do samtools index $file
done

date
