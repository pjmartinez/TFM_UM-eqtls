
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


module load hisat2/2.2.1  #load the module before uses it, very important.




for file in *_trim.fastq.gz
do name=$(basename $file _trim.fastq.gz) 

hisat2 -p 8 --dta -x /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/hisat_index/pinotnoir -U ${name}_trim.fastq.gz | \
        samtools view -S -h -u - | \
        samtools sort -T ${name} - > ./"$dir"/${name}.bam

echo "===================hisat align started at `date` for ============" $name
