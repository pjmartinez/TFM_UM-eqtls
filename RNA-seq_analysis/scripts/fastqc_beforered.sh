#!/bin/bash

#
#SBATCH -J fastqc_trimfiles_redvar
#SBATCH -p generic       # Partición (cola)
##SBATCH -N 1                # Numero de nodos
##SBATCH -n 8                 # Numero de cores(CPUs)
#SBATCH -w trueno188
#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 5-02:00     # Duración (D-HH:MM)
#SBATCH --output=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.out  #STDOUT
#SBATCH --error=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH --mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección

#export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch


cd /home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis
DIR=/home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/original_files

echo `hostname`

#################################################################
# FASTQC of raw reads 
#################################################################
dir="before"
if [ ! -d "$dir" ]; then
        mkdir -p $dir
fi


for file in ${DIR}/*.fastq.gz
do name=$(basename $file _Vitis_vinifera_RNA-Seq.fastq.gz)

/home/cebas/pmartinez/FastQC/fastqc --outdir ./"$dir"/ ${DIR}/${name}_Vitis_vinifera_RNA-Seq.fastq.gz -t 8




echo "===================fastqc original red varities at `date` for ============" $name

done
