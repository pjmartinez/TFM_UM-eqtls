
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




curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631879/SRR1631879.fastq.gz -o SRR1631879_GSM1532832_Vermentino_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631883/SRR1631883.fastq.gz -o SRR1631883_GSM1532836_Vermentino_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631880/SRR1631880.fastq.gz -o SRR1631880_GSM1532833_Vermentino_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631884/SRR1631884.fastq.gz -o SRR1631884_GSM1532837_Vermentino_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631881/SRR1631881.fastq.gz -o SRR1631881_GSM1532834_Vermentino_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631882/SRR1631882.fastq.gz -o SRR1631882_GSM1532835_Vermentino_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631886/SRR1631886.fastq.gz -o SRR1631886_GSM1532839_Vermentino_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631888/SRR1631888.fastq.gz -o SRR1631888_GSM1532841_Vermentino_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631889/SRR1631889.fastq.gz -o SRR1631889_GSM1532842_Vermentino_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631890/SRR1631890.fastq.gz -o SRR1631890_GSM1532843_Vermentino_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631885/SRR1631885.fastq.gz -o SRR1631885_GSM1532838_Vermentino_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631887/SRR1631887.fastq.gz -o SRR1631887_GSM1532840_Vermentino_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631892/SRR1631892.fastq.gz -o SRR1631892_GSM1532845_Garganega_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631894/SRR1631894.fastq.gz -o SRR1631894_GSM1532847_Garganega_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631891/SRR1631891.fastq.gz -o SRR1631891_GSM1532844_Garganega_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631893/SRR1631893.fastq.gz -o SRR1631893_GSM1532846_Garganega_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631895/SRR1631895.fastq.gz -o SRR1631895_GSM1532848_Garganega_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631896/SRR1631896.fastq.gz -o SRR1631896_GSM1532849_Garganega_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631897/SRR1631897.fastq.gz -o SRR1631897_GSM1532850_Garganega_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631898/SRR1631898.fastq.gz -o SRR1631898_GSM1532851_Garganega_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631899/SRR1631899.fastq.gz -o SRR1631899_GSM1532852_Garganega_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631900/SRR1631900.fastq.gz -o SRR1631900_GSM1532853_Garganega_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631902/SRR1631902.fastq.gz -o SRR1631902_GSM1532855_Garganega_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631901/SRR1631901.fastq.gz -o SRR1631901_GSM1532854_Garganega_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631903/SRR1631903.fastq.gz -o SRR1631903_GSM1532856_Glera_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631905/SRR1631905.fastq.gz -o SRR1631905_GSM1532858_Glera_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631904/SRR1631904.fastq.gz -o SRR1631904_GSM1532857_Glera_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631906/SRR1631906.fastq.gz -o SRR1631906_GSM1532859_Glera_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631908/SRR1631908.fastq.gz -o SRR1631908_GSM1532861_Glera_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631907/SRR1631907.fastq.gz -o SRR1631907_GSM1532860_Glera_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631910/SRR1631910.fastq.gz -o SRR1631910_GSM1532863_Glera_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631909/SRR1631909.fastq.gz -o SRR1631909_GSM1532862_Glera_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631911/SRR1631911.fastq.gz -o SRR1631911_GSM1532864_Glera_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631912/SRR1631912.fastq.gz -o SRR1631912_GSM1532865_Glera_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631913/SRR1631913.fastq.gz -o SRR1631913_GSM1532866_Glera_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631914/SRR1631914.fastq.gz -o SRR1631914_GSM1532867_Glera_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631915/SRR1631915.fastq.gz -o SRR1631915_GSM1532868_Moscatobianco_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631916/SRR1631916.fastq.gz -o SRR1631916_GSM1532869_Moscatobianco_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631917/SRR1631917.fastq.gz -o SRR1631917_GSM1532870_Moscatobianco_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631918/SRR1631918.fastq.gz -o SRR1631918_GSM1532871_Moscatobianco_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631920/SRR1631920.fastq.gz -o SRR1631920_GSM1532873_Moscatobianco_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631919/SRR1631919.fastq.gz -o SRR1631919_GSM1532872_Moscatobianco_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631923/SRR1631923.fastq.gz -o SRR1631923_GSM1532876_Moscatobianco_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631921/SRR1631921.fastq.gz -o SRR1631921_GSM1532874_Moscatobianco_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631922/SRR1631922.fastq.gz -o SRR1631922_GSM1532875_Moscatobianco_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631924/SRR1631924.fastq.gz -o SRR1631924_GSM1532877_Moscatobianco_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631926/SRR1631926.fastq.gz -o SRR1631926_GSM1532879_Moscatobianco_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631925/SRR1631925.fastq.gz -o SRR1631925_GSM1532878_Moscatobianco_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631927/SRR1631927.fastq.gz -o SRR1631927_GSM1532880_Passerina_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/009/SRR1631929/SRR1631929.fastq.gz -o SRR1631929_GSM1532882_Passerina_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631928/SRR1631928.fastq.gz -o SRR1631928_GSM1532881_Passerina_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/000/SRR1631930/SRR1631930.fastq.gz -o SRR1631930_GSM1532883_Passerina_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/001/SRR1631931/SRR1631931.fastq.gz -o SRR1631931_GSM1532884_Passerina_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/002/SRR1631932/SRR1631932.fastq.gz -o SRR1631932_GSM1532885_Passerina_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/003/SRR1631933/SRR1631933.fastq.gz -o SRR1631933_GSM1532886_Passerina_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/004/SRR1631934/SRR1631934.fastq.gz -o SRR1631934_GSM1532887_Passerina_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/005/SRR1631935/SRR1631935.fastq.gz -o SRR1631935_GSM1532888_Passerina_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/007/SRR1631937/SRR1631937.fastq.gz -o SRR1631937_GSM1532890_Passerina_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/006/SRR1631936/SRR1631936.fastq.gz -o SRR1631936_GSM1532889_Passerina_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR163/008/SRR1631938/SRR1631938.fastq.gz -o SRR1631938_GSM1532891_Passerina_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
