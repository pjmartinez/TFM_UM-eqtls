#!/bin/bash 
#SBATCH --job-name=deduplication
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=200G
#SBATCH --qos=himem
#SBATCH --partition=himem
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


module load bcftools


#INDIR=../align_stepwise
# make output directory if it doesn't exist. 
INDIR=../variants_bcftools
#mkdir -p $OUTDIR


# make a list of bam files
#ls $INDIR/*_mkdup.bam >$INDIR/list.bam

#GEN=/home/CAM/pmartinez/DATOSTFM/TFM_work/rawdata/PinotNoir.fa 

bcftools call -m -v -Oz -o $INDIR/vitis.vcf.gz $INDIR/vitis.pileup

date
