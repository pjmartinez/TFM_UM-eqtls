#!/bin/bash

#
#SBATCH -J HTcounts_red
#SBATCH -p generic       # Partición (cola)
##SBATCH -N 1                # Numero de nodos
##SBATCH -n 4                 # Numero de cores(CPUs)
#SBATCH -w trueno186
#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 5-02:00     # Duración (D-HH:MM)
#SBATCH --output=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.out  #STDOUT
#SBATCH --error=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH --mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección

#export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

############################
##fastqc and trimming files
############################


hostname
date

#Quality check of Reads
#DIR=/scratch-global/pmartinez/genotypes
#DIR=/home/cebas/pmartinez/secuencias/TFM_vitis/genotypes


DIR=/home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/trimmed_files/aligments_red


for file in ${DIR}/*.bam
do name=$(basename $file .bam)


echo "==========================counts started at `date` for ================" $name


htseq-count  -s no -r pos -f bam ${DIR}/${name}.bam /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/Vitis_vinifera_gene_annotation_on_V2_20_myversion.gff3  > ${DIR}/${name}.counts






done

