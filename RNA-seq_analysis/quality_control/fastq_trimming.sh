#!/bin/bash

#SBATCH --job-name=JOBNAME #Gives a user specified name to the job.
#SBATCH -n 1 #Task count
#SBATCH -N 1 #Node count
#SBATCH -c 1 #CPUs/cores per task
#SBATCH --mem=1G #job memory request per node, usually an integer followed by a prefix for the unit (e. g. --mem=1G for 1 GB)
#SBATCH --partition=general # Run the job in the specified partition/queue depend of your server.
#SBATCH --qos= general #Defines the quality-of-service to be used for the job.
#SBATCH --mail-type=ALL #Defines when a mail message about the job will be sent to the user. See the man page for details.
#SBATCH --mail-user=youremail
#SBATCH -o %x_%j.out #Specifies the file name to be used for stdout.
#SBATCH -e %x_%j.err #Specifies the file name to be used for stderr.


DIR=path to raw_data
DIR2= path to adapters

for file in ${DIR}/*.fastq.gz 
do name=$(basename $file _Vitis_vinifera_RNA-Seq.fastq.gz)

java -jar  /home/cebas/pmartinez/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        ${DIR}/${name}_Vitis_vinifera_RNA-Seq.fastq.gz  \
        ${DIR}/${name}_trim.fastq.gz \
        ILLUMINACLIP:${DIR2}/TruSeq2-SE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:45

#samtools sort ${DIR}/${name}.bam -o ${DIR}/${name}.sorted.bam
 
echo "===================trimming with trimomatic at `date` for ============" $name

done
