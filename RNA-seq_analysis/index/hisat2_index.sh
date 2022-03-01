
#!/bin/bash

#
#SBATCH -J hisat_index
#SBATCH -p generic       # Partición (cola)
##SBATCH -N 1                # Numero de nodos
##SBATCH -n 16                 # Numero de cores(CPUs)
#SBATCH -w trueno186
#SBATCH --mem-per-cpu=5000    # Bloque de memoria para todos los nodos
#SBATCH -t 5-02:00     # Duración (D-HH:MM)
#SBATCH --output=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.out  #STDOUT
#SBATCH --error=/home/cebas/pmartinez/secuencias/TFM_vitis/temp/slurm-%j.err   #STDERR
#SBATCH --mail-type=END,FAIL      # Notificación cuando el trabajo termina o falla
#SBATCH --mail-user=pjmartinez@cebas.csic.es # Enviar correo a la dirección


module load hisat2/2.2.1

hisat2-build -p 16 /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/all_chromosomes/PinotNoir.fa /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/hisat_index/pinotnoir
