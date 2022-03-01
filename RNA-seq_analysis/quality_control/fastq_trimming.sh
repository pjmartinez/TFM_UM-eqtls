#!/bin/bash

#
#SBATCH -J fastq_trimming
#SBATCH -p generic       # Partición (cola)
##SBATCH -N 1                # Numero de nodos
##SBATCH -n 1                 # Numero de cores(CPUs) 
#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 5-02:00     # Duración (D-HH:MM)
#SBATCH --output=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.out  #STDOUT
#SBATCH --error=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH --mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección



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
