#!/bin/bash 
#SBATCH --job-name=bwa_index
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=40G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=pjmartinezgarcia1@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date

module load bwa
module load samtools/1.10

cd ../genome
#GEN=/home/CAM/pmartinez/DATOSTFM/TFM_work/rawdata/GCF_000003745.3_12X_genomic.fa.gz

GEN=/home/CAM/pmartinez/DATOSTFM/TFM_work/rawdata/PinotNoir.fa

bwa index -p pinotnoir $GEN
#samtools faidx $GEN
