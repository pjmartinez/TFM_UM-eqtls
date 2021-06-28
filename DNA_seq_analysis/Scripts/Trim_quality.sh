#!/bin/bash 
#SBATCH --job-name=quality_control
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --qos=general
#SBATCH --partition=general
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o /home/CAM/pmartinez/DATOSTFM/TFM_work/scripts/temp/%x_%j.out
#SBATCH -e /home/CAM/pmartinez/DATOSTFM/TFM_work/scripts/temp/%x_%j.err

hostname
date

mkdir ../trimmed
#mkdir ../trimmed

module load sickle
module load fastqc


#Quality check of Reads

DIR=/home/CAM/pmartinez/DATOSTFM/TFM_work/rawdata

for file in ${DIR}/*_1.fastq.gz
do name=$(basename $file _1.fastq.gz)

#fastqc -t 4 -o ../fastqc ${DIR}/${name}_1.fastq.gz  ${DIR}/${name}_2.fastq.gz


#Trimming of reads
echo "========================================== Triming starting at `date` for  ==========================" $name

sickle pe \
-t sanger \
-f ${DIR}/${name}_1.fastq.gz \
-r ${DIR}/${name}_2.fastq.gz \
-o ../trimmed/${name}_trimmed_1.fastq.gz \
-p ../trimmed/${name}_trimmed_2.fastq.gz \
-l 45 -q 25 \
-s ../trimmed/singles_${name}.fastq.gz


echo "========================================== Triming finished at `date` for  ==========================" $name
echo "========================================== fastqc starting at `date` for  ==========================" $name


#Quality check of Reads post trimming
fastqc -t 4 -o ../fastqc ../trimmed/${name}_trimmed_1.fastq.gz ../trimmed/${name}_trimmed_2.fastq.gz

echo "===========================================fastqc finished at `date` for  ==========================="$name


done
