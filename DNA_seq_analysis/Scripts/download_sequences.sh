#!/bin/bash 
#SBATCH --job-name=align
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=40G
#SBATCH --qos=special
#SBATCH --partition=special_himem
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

hostname
date


cd ../rawdata

#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156354/SRR6156354_1.fastq.gz -o SRR6156354_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156354/SRR6156354_2.fastq.gz -o SRR6156354_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/006/SRR6156356/SRR6156356_1.fastq.gz -o SRR6156356_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/006/SRR6156356/SRR6156356_2.fastq.gz -o SRR6156356_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/004/SRR5506714/SRR5506714_1.fastq.gz -o SangioveseC_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/004/SRR5506714/SRR5506714_2.fastq.gz -o SangioveseC_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/005/SRR5506715/SRR5506715_1.fastq.gz -o SangioveseB_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/005/SRR5506715/SRR5506715_2.fastq.gz -o SangioveseB_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/001/SRR5506711/SRR5506711_1.fastq.gz -o SangioveseF_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/001/SRR5506711/SRR5506711_2.fastq.gz -o SangioveseF_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/003/SRR5506713/SRR5506713_1.fastq.gz -o SangioveseD_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/003/SRR5506713/SRR5506713_2.fastq.gz -o SangioveseD_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/002/SRR5506712/SRR5506712_1.fastq.gz -o SangioveseE_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/002/SRR5506712/SRR5506712_2.fastq.gz -o SangioveseE_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/006/SRR5506716/SRR5506716_1.fastq.gz -o SangioveseA_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/006/SRR5506716/SRR5506716_2.fastq.gz -o SangioveseA_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/007/SRR5506717/SRR5506717_1.fastq.gz -o CabernetSauvignon_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/007/SRR5506717/SRR5506717_2.fastq.gz -o CabernetSauvignon_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156334/SRR6156334_1.fastq.gz -o SRR6156334_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156334/SRR6156334_2.fastq.gz -o SRR6156334_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/007/SRR6156337/SRR6156337_1.fastq.gz -o SRR6156337_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/007/SRR6156337/SRR6156337_2.fastq.gz -o SRR6156337_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR580/008/SRR5803838/SRR5803838_1.fastq.gz -o Vermentino_resequencing_1.fastq.gz
#curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR580/008/SRR5803838/SRR5803838_2.fastq.gz -o Vermentino_resequencing_2.fastq.gz
#


curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/000/SRR6156340/SRR6156340_1.fastq.gz -o SRR6156340_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/000/SRR6156340/SRR6156340_2.fastq.gz -o SRR6156340_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/003/SRR6156343/SRR6156343_1.fastq.gz -o SRR6156343_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/003/SRR6156343/SRR6156343_2.fastq.gz -o SRR6156343_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/002/SRR6156342/SRR6156342_1.fastq.gz -o SRR6156342_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/002/SRR6156342/SRR6156342_2.fastq.gz -o SRR6156342_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/001/SRR6156341/SRR6156341_1.fastq.gz -o SRR6156341_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/001/SRR6156341/SRR6156341_2.fastq.gz -o SRR6156341_DNA-Seq_of_Vitis_vinifera_2.fastq.gz




