






for file in *_trim.fastq.gz
do name=$(basename $file _trim.fastq.gz) 

hisat2 -p 8 --dta -x /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/hisat_index/pinotnoir -U ${name}_trim.fastq.gz | \
        samtools view -S -h -u - | \
        samtools sort -T ${name} - > ./"$dir"/${name}.bam

echo "===================hisat align started at `date` for ============" $name
