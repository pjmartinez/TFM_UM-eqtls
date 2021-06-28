#!/bin/bash 
#SBATCH --job-name=bam_sort
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH --qos=general
#SBATCH --partition=amd
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch


cd ../align

for file in *.bam
do name=$(basename $file .bam) 

echo "===================sort bam  `date` for ============" $name

java -jar $PICARD SortSam \
        INPUT=${name}.bam \
        OUTPUT=../align_stepwise/${name}_sort.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True

done
