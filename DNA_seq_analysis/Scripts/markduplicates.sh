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
DIR2=/home/CAM/pmartinez/DATOSTFM/TFM_work/align_stepwise

module load picard/2.9.2
export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch


for file in ${DIR2}/*_sort.bam
do name=$(basename $file _sort.bam)

java -jar $PICARD MarkDuplicates \
        INPUT=${DIR2}/${name}_sort.bam \
        OUTPUT=${DIR2}/${name}_mkdup.bam \
        REMOVE_DUPLICATES=Ture \
        METRICS_FILE=${DIR2}/${name}_mkdup_metrics.txt \
        CREATE_INDEX=True
done
