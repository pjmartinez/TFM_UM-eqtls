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

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631819/SRR1631819.fastq.gz -o SRR1631819_GSM1532772_Sangiovese_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631822/SRR1631822.fastq.gz -o SRR1631822_GSM1532775_Sangiovese_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631821/SRR1631821.fastq.gz -o SRR1631821_GSM1532774_Sangiovese_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631820/SRR1631820.fastq.gz -o SRR1631820_GSM1532773_Sangiovese_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631824/SRR1631824.fastq.gz -o SRR1631824_GSM1532777_Sangiovese_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631823/SRR1631823.fastq.gz -o SRR1631823_GSM1532776_Sangiovese_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631826/SRR1631826.fastq.gz -o SRR1631826_GSM1532779_Sangiovese_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631827/SRR1631827.fastq.gz -o SRR1631827_GSM1532780_Sangiovese_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631828/SRR1631828.fastq.gz -o SRR1631828_GSM1532781_Sangiovese_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631825/SRR1631825.fastq.gz -o SRR1631825_GSM1532778_Sangiovese_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631829/SRR1631829.fastq.gz -o SRR1631829_GSM1532782_Sangiovese_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631830/SRR1631830.fastq.gz -o SRR1631830_GSM1532783_Sangiovese_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631831/SRR1631831.fastq.gz -o SRR1631831_GSM1532784_Barbera_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631833/SRR1631833.fastq.gz -o SRR1631833_GSM1532786_Barbera_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631834/SRR1631834.fastq.gz -o SRR1631834_GSM1532787_Barbera_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631832/SRR1631832.fastq.gz -o SRR1631832_GSM1532785_Barbera_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631835/SRR1631835.fastq.gz -o SRR1631835_GSM1532788_Barbera_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631836/SRR1631836.fastq.gz -o SRR1631836_GSM1532789_Barbera_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631837/SRR1631837.fastq.gz -o SRR1631837_GSM1532790_Barbera_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631838/SRR1631838.fastq.gz -o SRR1631838_GSM1532791_Barbera_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631839/SRR1631839.fastq.gz -o SRR1631839_GSM1532792_Barbera_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631840/SRR1631840.fastq.gz -o SRR1631840_GSM1532793_Barbera_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631841/SRR1631841.fastq.gz -o SRR1631841_GSM1532794_Barbera_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631842/SRR1631842.fastq.gz -o SRR1631842_GSM1532795_Barbera_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631843/SRR1631843.fastq.gz -o SRR1631843_GSM1532796_Negroamaro_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631844/SRR1631844.fastq.gz -o SRR1631844_GSM1532797_Negroamaro_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631845/SRR1631845.fastq.gz -o SRR1631845_GSM1532798_Negroamaro_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631846/SRR1631846.fastq.gz -o SRR1631846_GSM1532799_Negroamaro_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631847/SRR1631847.fastq.gz -o SRR1631847_GSM1532800_Negroamaro_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631848/SRR1631848.fastq.gz -o SRR1631848_GSM1532801_Negroamaro_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631849/SRR1631849.fastq.gz -o SRR1631849_GSM1532802_Negroamaro_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631851/SRR1631851.fastq.gz -o SRR1631851_GSM1532804_Negroamaro_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631852/SRR1631852.fastq.gz -o SRR1631852_GSM1532805_Negroamaro_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631850/SRR1631850.fastq.gz -o SRR1631850_GSM1532803_Negroamaro_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631853/SRR1631853.fastq.gz -o SRR1631853_GSM1532806_Negroamaro_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631854/SRR1631854.fastq.gz -o SRR1631854_GSM1532807_Negroamaro_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631855/SRR1631855.fastq.gz -o SRR1631855_GSM1532808_Refosco_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631856/SRR1631856.fastq.gz -o SRR1631856_GSM1532809_Refosco_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631857/SRR1631857.fastq.gz -o SRR1631857_GSM1532810_Refosco_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631858/SRR1631858.fastq.gz -o SRR1631858_GSM1532811_Refosco_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631859/SRR1631859.fastq.gz -o SRR1631859_GSM1532812_Refosco_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631860/SRR1631860.fastq.gz -o SRR1631860_GSM1532813_Refosco_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631861/SRR1631861.fastq.gz -o SRR1631861_GSM1532814_Refosco_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631862/SRR1631862.fastq.gz -o SRR1631862_GSM1532815_Refosco_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631863/SRR1631863.fastq.gz -o SRR1631863_GSM1532816_Refosco_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631864/SRR1631864.fastq.gz -o SRR1631864_GSM1532817_Refosco_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631865/SRR1631865.fastq.gz -o SRR1631865_GSM1532818_Refosco_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631866/SRR1631866.fastq.gz -o SRR1631866_GSM1532819_Refosco_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631867/SRR1631867.fastq.gz -o SRR1631867_GSM1532820_Primitivo_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631868/SRR1631868.fastq.gz -o SRR1631868_GSM1532821_Primitivo_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631870/SRR1631870.fastq.gz -o SRR1631870_GSM1532823_Primitivo_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631869/SRR1631869.fastq.gz -o SRR1631869_GSM1532822_Primitivo_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631871/SRR1631871.fastq.gz -o SRR1631871_GSM1532824_Primitivo_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631872/SRR1631872.fastq.gz -o SRR1631872_GSM1532825_Primitivo_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631873/SRR1631873.fastq.gz -o SRR1631873_GSM1532826_Primitivo_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631874/SRR1631874.fastq.gz -o SRR1631874_GSM1532827_Primitivo_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631877/SRR1631877.fastq.gz -o SRR1631877_GSM1532830_Primitivo_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631878/SRR1631878.fastq.gz -o SRR1631878_GSM1532831_Primitivo_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631876/SRR1631876.fastq.gz -o SRR1631876_GSM1532829_Primitivo_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631875/SRR1631875.fastq.gz -o SRR1631875_GSM1532828_Primitivo_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
