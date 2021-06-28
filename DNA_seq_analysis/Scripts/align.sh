#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=40G
#SBATCH --qos=special
#SBATCH --partition=special_himem
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load bwa/0.7.17
module load samtools/1.7

#mkdir ../align

cd /home/CAM/pmartinez/DATOSTFM/TFM_work/trimmed

GEN=/home/CAM/pmartinez/DATOSTFM/TFM_work/genome/pinotnoir



for file in *_1.fastq.gz 
do name=$(basename $file _trimmed_1.fastq.gz) 

echo "===================bwa -mem started at `date` for ============" $name
bwa mem -t 12  $GEN  ${name}_trimmed_1.fastq.gz ${name}_trimmed_2.fastq.gz  -o ../align/${name}.sam

echo "===================bwa -mem finished at `date` for ============" $name
#SAM to BAM CONVERSION

echo "===================sam started at `date` for ============" $name

samtools view -bhS ../align/${name}.sam > ../align/${name}.bam

echo "===================sam finished at `date` for ============" $name

done
