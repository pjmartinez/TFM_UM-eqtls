#!/bin/bash
#
#SBATCH -p generic       # Partición (cola)
#SBATCH -N 1                # Numero de nodos
#SBATCH -n 1                 # Numero de cores(CPUs)

#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 0-02:00     # Duración (D-HH:MM)
#SBATCH -o slurm.%N, %j.out  #STDOUT
#SBATCH -e slurm.%N, %j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH -- mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección


cd /home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/trimmed_files/aligments_red

for file in *.bam
do

/home/cebas/pmartinez/samtools-1.10/samtools index $file

done

