#!/bin/bash 
#SBATCH --job-name=samtools
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load bcftools
module load samtools


DIR2=/home/CAM/pmartinez/DATOSTFM/TFM_work/align_stepwise



for file in ${DIR2}/*_mkdup.bam
do name=$(basename $file _mkdup.bam)
echo $file


# make a list of bam files

GEN=/home/CAM/pmartinez/DATOSTFM/TFM_work/rawdata/PinotNoir.fa 

samtools stats -r $GEN ${DIR2}/${name}_mkdup.bam > ${DIR2}/${name}_samstat.txt

date
done
