#!/bin/bash

#
#SBATCH -J alig_hisat_red
#SBATCH -p generic       # Partición (cola)
##SBATCH -N 1                # Numero de nodos
##SBATCH -n 8                 # Numero de cores(CPUs)
#SBATCH -w trueno186
#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 5-02:00     # Duración (D-HH:MM)
#SBATCH --output=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.out  #STDOUT
#SBATCH --error=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH --mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección

#export _JAVA_OPTIONS=-Djava.io.tmpdir=/scratch

#mkdir ../align_stepwise
#DIR1=/home/cebas/pmartinez/sam_files_pomogranate
#SCRATCHDIR=/scratch-global/$USER/$SLURM_JOBID
#DIR1=$SCRATCHDIR
#DIR2=/home/cebas/pmartinez/Sara_pomogranates/align_stepwise

cd /home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/trimmed_files


echo `hostname`


############################
##Hisat for white vitis
############################

dir="aligments_red"

if [ ! -d "$dir" ]; then
        mkdir -p $dir
fi



for file in *_trim.fastq.gz
do name=$(basename $file _trim.fastq.gz) 

/home/cebas/pmartinez/hisat2-2.1.0/hisat2 -p 8 --dta -x /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/hisat_index/pinotnoir -U ${name}_trim.fastq.gz | \
        samtools view -S -h -u - | \
        samtools sort -T ${name} - > ./"$dir"/${name}.bam

echo "===================hisat align started at `date` for ============" $name



done
