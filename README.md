
# Bioinformatic approach for the discovery of cis-eQTL signals during fruit ripening of a woody species as grape (*Vitis vinifera* L.).


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a [Slurm](https://slurm.schedmd.com/sbatch.html) scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. More information about [Slurm](https://slurm.schedmd.com/sbatch.html) can be found in the [hpc wiki](https://hpc-wiki.info/hpc/SLURM). Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
In summary, the repository includs the different inputs, scripts, outputs files (the supplemental files) obtained during the analysis performance in the biorxv paper: https://www.biorxiv.org/content/10.1101/2021.07.06.450811v1





![Screenshot](/Figures/generalpipeline.png)

## Cloning the tutorial repository

To work through this tutorial, copy it to your home directory using the git clone command:
```
git clone < git-repository-path.git > 

```

Contents
 1. [Overview](#overview)
 2. [RNA Seq](#rna-seq-analysis)
 3. [Variant Calling](#Variant-Calling)
 4. [eQTL analysis](#eQTL-analysis)
 5. [GO ontology annotation with biomaRt](#GO-biomart)

# Overview
This tutorial will teach you how to use open source quality control, RNA Seq, Variant Calling, eQTL tools to complete a cis-eQTL analysis which is possible when you  have generated the specific datasets. Moving through the tutorial, you will take expression and genotypic data from a woody species as grape and perform a eQTL analysis via Matrix eQTL to characterize the gene expression levels during ripening Vitis vinifera L. fruit.

# RNA Seq Analysis
A total of 400mg of RNA was extracted from berry pericarp tissue (entire berries without seeds), a detailed description about RNA extraction and library preparation and sequencing of the 120 samples (10 varieties at four stages, in total 40 triplicate samples) can be also found in [Massonnet, M. et al. 2017](https://pubmed.ncbi.nlm.nih.gov/28652263/). 




## 2.1. Accessing the Data using SRA-Toolkit
Before we can get started, we need to get the data we're going to analyze. This dataset has been deposited in the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) at NCBI, a comprehensive collection of sequenced genetic data submitted by researchers. The beauty of the SRA is the ease with which genetic data becomes accessible to any scientist with an internet connection. Sets of sequences (usually all the sequences from a given sample within an experiment) in the SRA have a unique identifier. The set may be downloaded using a software module called the `sra-explorer`. There are several possibilities to download the files in the `sra-explorer`, which I invite you to investigate for yourself at [here](https://sra-explorer.info/).

An overview of the project for the RNA data can be viewed [here](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA265039) and [here](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA265040).

We will download data from all growtn stages (Pea-sized, Prior veraison, end of veraison and Harvest). The SRA accessions are as follows:

```
Run	growth_stage
SRR1631819	Pea-sized berries
SRR1631820	Pea-sized berries
SRR1631821	Pea-sized berries
SRR1631822	Berries beginning to touch
SRR1631823	Berries beginning to touch
SRR1631824	Berries beginning to touch
SRR1631825	Soft berries
SRR1631826	Soft berries
SRR1631827	Soft berries
SRR1631828	Berries ripe for harvest
SRR1631829	Berries ripe for harvest
SRR1631830	Berries ripe for harvest
SRR1631831	Pea-sized berries
SRR1631832	Pea-sized berries
SRR1631833	Pea-sized berries
SRR1631834	Berries beginning to touch
SRR1631835	Berries beginning to touch
SRR1631836	Berries beginning to touch
.
.
.
.
.


```

We have provided a script to download the data from the the SRA data using this script. It contains two sections.

The SLURM header:
```
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

The download commands:

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
.
.
.

```

The full script for slurm scheduler for red and white cultivars can be found in the [raw_data](https://github.com/pjmartinez/TFM_UM-eqtls/tree/main/rawdata) folder. Before running it, add your own e-mail address to the --mail-user option (or delete the line entirely if you don't want an e-mail notification when the job completes).

When you're ready, you can execute the script by entering sbatch [rna_red_download.sh](https://github.com/pjmartinez/TFM_UM-eqtls/tree/main/rawdata/rna_red_download.sh) or sbatch [rna_white_download.sh](https://github.com/pjmartinez/TFM_UM-eqtls/tree/main/rawdata/rna_white_download.sh) in the terminal. This submits the jobs to the SLURM scheduler.

Once the job is completed the folder structure will look like this:

raw_data/
```
SRR1631819_GSM1532772_Sangiovese_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631822_GSM1532775_Sangiovese_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631821_GSM1532774_Sangiovese_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631820_GSM1532773_Sangiovese_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631824_GSM1532777_Sangiovese_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631823_GSM1532776_Sangiovese_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631826_GSM1532779_Sangiovese_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631827_GSM1532780_Sangiovese_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631828_GSM1532781_Sangiovese_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631825_GSM1532778_Sangiovese_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631829_GSM1532782_Sangiovese_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631830_GSM1532783_Sangiovese_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631831_GSM1532784_Barbera_Pea_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631833_GSM1532786_Barbera_Pea_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631834_GSM1532787_Barbera_Touch_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631832_GSM1532785_Barbera_Pea_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631835_GSM1532788_Barbera_Touch_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631836_GSM1532789_Barbera_Touch_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631837_GSM1532790_Barbera_Soft_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631838_GSM1532791_Barbera_Soft_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631839_GSM1532792_Barbera_Soft_3_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631840_GSM1532793_Barbera_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631841_GSM1532794_Barbera_Harv_2_Vitis_vinifera_RNA-Seq.fastq.gz
SRR1631842_GSM1532795_Barbera_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz
...
```
In both scripts are run a total of 120 files should be there. The sequence files are in fastq format and compressed using gzip (indicated by the .gz). It's good practice to keep sequence files compressed. Most bioinformatics programs can read them directly without needing to decompress them first, and it doesn't get in the way of inspecting them either. The .out and .err files are output produced by SLURM that you can use to troubleshoot if things go wrong. Lets have a look at at the contents of one of the fastq files:

zcat SRR1631842_GSM1532795_Barbera_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz | head -n 12
```
@SRR1631844.1 HWI-1KL152:93:D0YY9ACXX:2:1101:1175:2042/1
NAGGGCCCTGCTGAAGGAATTTTGATTTCTGGTTACTCACTTGCAATTGATGAATCAAGAATAACTGGATATAGCAAGAATGTTCATAATGAAATA
+
################################################################################################
@SRR1631844.2 HWI-1KL152:93:D0YY9ACXX:2:1101:1136:2046/1
NCCAGATTANNCNNNNNNNCATTCACGTTAGTTTCTGANGNANTANCAGGCCCTNTATTCTCATCTGTTTNNCTCTTTANGCGACANAGNNNANNG
+
#0;9@@?=?##2#######32@@@<?@1=????=?#############################################################
@SRR1631844.3 HWI-1KL152:93:D0YY9ACXX:2:1101:1175:2084/1
CTTTGACAGTCTTTGCTTCCTTGATGGCAGCACGAATCTCATCATAGCCAGTATTTCAATTCTTCACCCAGATAACATGACAACCAAGAGCCTCAA
+
;;1((@??<.4=:@=)2;))86);)9))=((3:8?:>?=?=>>>8>))07>7>????#######################################

```

Each sequence record has four lines. The first is the sequence name, beginning with @. The second is the nucleotide sequence. The third is a comment line, beginning with +, and which here only contains the sequence name again (it is often empty). The fourth are phred-scaled base quality scores, encoded by ASCII characters. Follow the links to learn more, but in short, the quality scores give the probability a called base is incorrect.



## 2.2. Quality Control Using Trimmomatic
`Trimmomatic` is commonly used to trim low quality and adapter contaminated sequences.

Our usage looks like this for a single sample:
```

module load Trimmomatic/0.39

DIR=path to raw_data
DIR2= path to adapters

for file in ${DIR}/*.fastq.gz #loop for search each gz file
do name=$(basename $file _Vitis_vinifera_RNA-Seq.fastq.gz) # to obtain a name for each sample used

java -jar  trimmomatic-0.39.jar SE \ #basic command for trimmomatic
        ${DIR}/${name}_Vitis_vinifera_RNA-Seq.fastq.gz  \
        ${DIR}/${name}_trim.fastq.gz \
        ILLUMINACLIP:${DIR2}/TruSeq2-SE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:45

```
 
We call SE for single-end mode and we specify the input and output file names. The ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 command searches for adapter sequence, so we provide a fasta file containing the adapters used in the library preparation, and the numbers control the parameters of adapter matching (see the manual for more details). SLIDINGWINDOW:4:20 scans through the read, cutting the read when the average base quality in a 4 base window drops below 20. We linked to an explanation of phred-scaled quality scores above, but for reference, scores of 10 and 20 correspond to base call error probabilities of 0.1 and 0.01, respectively. MINLEN:45 causes reads to be dropped if they have been trimmed to less than 45bp.Here is a useful paper on setting trimming parameters for RNA-seq.

The full scripts for slurm scheduler is calling which can be found in the [quality_control/](https://github.com/pjmartinez/TFM_UM-eqtls/tree/main/RNA-seq_analysis/quality_control) folder. Navigate there and run the script by enteriing sbatch [fastq_trimming.sh](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/RNA-seq_analysis/quality_control/fastq_trimming.sh) on the command-line.

Following the trimmomatic run, the resulting file structure will look as follows:

quality_control/

```
Sangiovese_Pea_1_trim.fastq.gz
Sangiovese_Touch_1_trim.fastq.gz
Sangiovese_Pea_3_trim.fastq.gz
Sangiovese_Pea_2_trim.fastq.gz
Sangiovese_Touch_3_trim.fastq.gz
Sangiovese_Touch_2_trim.fastq.gz
Sangiovese_Soft_2_trim.fastq.gz
Sangiovese_Soft_3_trim.fastq.gz
Sangiovese_Harv_1_trim.fastq.gz
Sangiovese_Soft_1_trim.fastq.gz
Sangiovese_Harv_2_trim.fastq.gz
Sangiovese_Harv_3_trim.fastq.gz
Barbera_Pea_1_trim.fastq.gz
Barbera_Pea_3_trim.fastq.gz
Barbera_Touch_1_trim.fastq.gz
Barbera_Pea_2_trim.fastq.gz
Barbera_Touch_2_trim.fastq.gz
Barbera_Touch_3_trim.fastq.gz
Barbera_Soft_1_trim.fastq.gz
Barbera_Soft_2_trim.fastq.gz
Barbera_Soft_3_trim.fastq.gz
Barbera_Harv_1_trim.fastq.gz
Barbera_Harv_2_trim.fastq.gz
Barbera_Harv_3_trim.fastq.gz
Negroamaro_Pea_1_trim.fastq.gz
Negroamaro_Pea_2_trim.fastq.gz
Negroamaro_Pea_3_trim.fastq.gz
Negroamaro_Touch_1_trim.fastq.gz
Negroamaro_Touch_2_trim.fastq.gz
Negroamaro_Touch_3_trim.fastq.gz
Negroamaro_Soft_1_trim.fastq.gz
Negroamaro_Soft_3_trim.fastq.gz
Negroamaro_Harv_1_trim.fastq.gz
Negroamaro_Soft_2_trim.fastq.gz
Negroamaro_Harv_2_trim.fastq.gz
Negroamaro_Harv_3_trim.fastq.gz
Refosco_Pea_1_trim.fastq.gz
Refosco_Pea_2_trim.fastq.gz
Refosco_Pea_3_trim.fastq.gz
Refosco_Touch_1_trim.fastq.gz
Refosco_Touch_2_trim.fastq.gz
Refosco_Touch_3_trim.fastq.gz
Refosco_Soft_1_trim.fastq.gz
Refosco_Soft_2_trim.fastq.gz
Refosco_Soft_3_trim.fastq.gz
Refosco_Harv_1_trim.fastq.gz
Refosco_Harv_2_trim.fastq.gz
Refosco_Harv_3_trim.fastq.gz
Primitivo_Pea_1_trim.fastq.gz
Primitivo_Pea_2_trim.fastq.gz
Primitivo_Touch_1_trim.fastq.gz
Primitivo_Pea_3_trim.fastq.gz
Primitivo_Touch_2_trim.fastq.gz
Primitivo_Touch_3_trim.fastq.gz
Primitivo_Soft_1_trim.fastq.gz
Primitivo_Soft_2_trim.fastq.gz
Primitivo_Harv_2_trim.fastq.gz
Primitivo_Harv_3_trim.fastq.gz
Primitivo_Harv_1_trim.fastq.gz
Primitivo_Soft_3_trim.fastq.gz
.
.
.


```


Examine the .out file generated during the run. Summaries of how many reads were retained for each file were written there. Here's one example:

```

TrimmomaticSE: Started with arguments:
 /home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/Barbera_Harv_1_Vitis_vinifera_RNA-Seq.fastq.gz /home/cebas/pmartinez/secuencias/TFM_vitis/RNA_seq_red_vitis/Barbera_Harv_1_trim.fastq.gz ILLUMINACLIP:/home/cebas/pmartinez/Trimmomatic-0.39/adapters/TruSeq2-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:45
Automatically using 1 threads
Using Long Clipping Sequence: 'AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG'
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
Using Long Clipping Sequence: 'AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG'
ILLUMINACLIP: Using 0 prefix pairs, 3 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Reads: 35848728 Surviving: 31761232 (88,60%) Dropped: 4087496 (11,40%)
TrimmomaticSE: Completed successfully
```

## 2.3. FASTQC Before and After Quality Control
It is helpful to see how the quality of the data has changed after using Trimmomatic. To do this, we will be using the command-line version of fastqc . This program creates visual reports of the average quality of our reads.

```

module load fastqc/0.11.7
cd path/to/RNAfolder
DIR=/path/to/rawdata/  # or any other path where raw data is located

echo `hostname`

#################################################################
# FASTQC of raw reads 
#################################################################
dir="before" #to generate a dir with data quality before trimming
if \[ ! -d "$dir" ]; then #if the dir no exist it make it now
        mkdir -p $dir
fi


for file in ${DIR}/*.fastq.gz

do name=$(basename $file _Vitis_vinifera_RNA-Seq.fastq.gz)

fastqc --outdir ./"$dir"/ ${DIR}/${name}_Vitis_vinifera_RNA-Seq.fastq.gz -t 8



echo "===================fastqc original grape varities at `date` for ============" $name

done
```

The same command can be run on the fastq files after the trimming using fastqc program, and the comand will look like this:

```
module load fastqc/0.11.7
cd path/to/RNAfolder
DIR=/path/to/qualitycontrol # or any other path where the trimmed files are located

echo `hostname`

#################################################################
# FASTQC of raw reads 
#################################################################
dir="after" to generate a dir with data quality after trimming
if \[ ! -d "$dir" ]; then
        mkdir -p $dir
fi


for file in ${DIR}/*.fastq.gz
do name=$(basename $file _trim.fastq.gz)

/home/cebas/pmartinez/FastQC/fastqc --outdir ./"$dir"/ ${DIR}/${name}_trim.fastq.gz -t 8




echo "===================fastqc trimmed grape varities at `date` for ============" $name

done

```

This will produce html files with the quality reports. The file structure inside the folder fastqc/ will look like this:

```         
           
fastqc/
├── after
│   ├── Garganega_Harv_1_trim_fastqc.html
│   ├── Garganega_Harv_1_trim_fastqc.zip
│   ├── Glera_Harv_1_trim_fastqc.html
│   ├── Glera_Harv_1_trim_fastqc.zip
│   ├── Passerina_Harv_1_trim_fastqc.html
│   └── Passerina_Harv_1_trim_fastqc.zip 
├── before
│   ├── Garganega_Harv_1_Vitis_vinifera_RNA-Seq_fastqc.html
│   ├── Garganega_Harv_1_Vitis_vinifera_RNA-Seq_fastqc.zip
│   ├── Glera_Soft_3_Vitis_vinifera_RNA-Seq_fastqc.html
│   ├── Glera_Soft_3_Vitis_vinifera_RNA-Seq_fastqc.zip 
│   ├── Passerina_Pea_2_Vitis_vinifera_RNA-Seq_fastqc.html
│   └── Passerina_Pea_2_Vitis_vinifera_RNA-Seq_fastqc.zip

```
To view the html files you need to download them to your laptop and open them in a web browser. You can get them [here](https://bk-genomica.cebas.csic.es:5001/sharing/fVxtvpfSz)

There are some basic statistics which are all pretty self-explanatory. You should notice that none of our sequence libraries fail the quality report! It would be concerning if we had even one because this report is from our trimmed sequence! The same thinking applies to our sequence length. Should the minimum of the sequence length be below 45, we would know that Trimmomatic had not run properly. 


## 2.4 Aligning Reads to a Genome using `HISAT2`

### Downloading the genome and building the Index:

HISAT2 is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. In order to map the reads to a reference genome we have to do a few things to prepare. First we must download the reference genome! We will download the reference genome of the cultivar [Pinot Noir PN40024] (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/745/GCF_000003745.3_12X/) from the ncbi database using the `wget` command.

```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/745/GCF_000003745.3_12X/GCF_000003745.3_12X_genomic.fna.gz
gunzip GCF_000003745.3_12X_genomic.fna.gz
```
Next, we need to create a genome *index*. What is an index and why is it helpful? Genome indexing is the same as indexing a tome, like an encyclopedia. It is much easier to locate information in the vastness of an encyclopedia when you consult the index, which is ordered in an easily navigable way with pointers to the information you seek within. Genome indexing similarly structures the information contained in a genome so that a read mapper can quickly find possible mapping locations.

We will use the `hisat2-build` module to make a HISAT index file for the genome. It will create a set of files with the suffix .ht2, these files together comprise the index. The command to generate the index looks like this:

```

module load hisat2/2.2.1
hisat2-build -p 16 path/to/genome nameoftheindex


```

The full script can be found in the **index** folder by the name [hisat2_index.sh](https://github.com/pjmartinez/TFM_UM-eqtls/tree/main/RNA-seq_analysis/index). Navigate there and submit it by entering sbatch hisat2_index.sh on the command-line.

After running the script, the following files will be generated as part of the index. To refer to the index for mapping the reads in the next step, you will use the file prefix, which in this case is: pinotnoir

```
index/
|-- pinotnoir.1.ht2
|-- pinotnoir.2.ht2
|-- pinotnoir.3.ht2
|-- pinotnoir.4.ht2
|-- pinotnoir.5.ht2
|-- pinotnoir.6.ht2
|-- pinotnoir.7.ht2
`-- pinotnoir.8.ht2
```

### Aligning the reads using HISAT2
Once we have created the index, the next step is to align the reads to the reference genome with `HISAT2`. By default `HISAT2` outputs the alignments in SAM format. We won't go over the format in detail in this tutorial, but should you actually need to look at the alignment files, it would be helpful to read over the [format specification](https://samtools.github.io/hts-specs/SAMv1.pdf) or have a look the [wikipedia page](https://en.wikipedia.org/wiki/SAM_(file_format)).

Raw `SAM` formatted files have two issues. First, they are uncompressed. Because of this they take up much more space than they need to and are slower to read and write. Second, `HISAT2` writes the alignments out in the same order it reads the sequences from the fastq file, but for downstream applications, they need to be sorted by **genome coordinate**.

To deal with these issues, we'll use a pipe to send the results from `HISAT2` to `samtools` to sort the sequences, convert them to binary format and compress them. The resulting files will be in BAM format. We can then use `samtools`to read or otherwise manipulate them. We use pipes rather than writing intermediate files because it is much more efficient computationally, and requires less cleanup of unneeded intermediate files.

Here's our code for aligning one sample:
```
module load hisat2/2.2.1
module load samtools/1.10
hisat2 -p 8 --dta -x /home/cebas/pmartinez/secuencias/TFM_vitis/PinorNoir_genome/hisat_index/pinotnoir -U ${name}_trim.fastq.gz | \
        samtools view -S -h -u - | \
        samtools sort -T ${name} - > ./"$dir"/${name}.bam

```
The `|` is the pipe. It tells linux to use the output of the command to the left of the pipe as the input for the command to the right. You can chain together many commands using pipes. `samtools` view converts the SAM file produced by `hisat2` to uncompressed BAM. `-S` indicates the input is SAM format. `-h` indicates the SAM header should be written to the output. `-u` indicates that uncompressed BAM format should be written (no need to compress until the end). `-` indicates samtools should take input from the pipe. `samtools sort` sorts and compressed the file. `-T` gives a temporary file prefix.

Because BAM files are large and we may want to access specific sections quickly, we need to index the bam files, just like we indexed the genome. Here we'll use `samtools` again. As an example:


```

module load samtools/1.10

for file in *.bam
do

samtools index $file

done
```

This creates a `.bam.bai` index file to accompany each BAM file.
You can acces to all `.bam` and `.bam.bai` files generated [here](https://bk-genomica.cebas.csic.es:5001/sharing/40ME2HGUc)

The full script for the slurm scheduler can be found in the **align/** directory by the name [align.sh](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/RNA-seq_analysis/align/aligh.sh). When you're ready, navigate there and execute it by entering `sbatch align.sh` on the command-line.

When HISAT2 finishes aligning all the reads, it will write a summary which will be captured by SLURM in the file ending .err.

An alignment summary for a single sample is shown below:
```

25664909 reads; of these:
  25664909 (100.00%) were unpaired; of these:
    1114878 (4.34%) aligned 0 times
    23209585 (90.43%) aligned exactly 1 time
    1340446 (5.22%) aligned >1 times
95.66% overall alignment rate

```

Let's have a look at the BAM file:
```
module load samtools/1.9
samtools view -H Barbera_Harv_1.bam | head

```

which will give:
```

@HD	VN:1.0	SO:coordinate
@SQ	SN:NC_012016.3	LN:18140952
@SQ	SN:NC_012017.3	LN:19818926
@SQ	SN:NC_012018.3	LN:22702307
@SQ	SN:NC_012019.3	LN:24396255
@SQ	SN:NC_012020.3	LN:30274277
@SQ	SN:NC_012021.3	LN:20304914
@SQ	SN:NC_012022.3	LN:22053297
@SQ	SN:NC_012023.3	LN:17126926
@SQ	SN:NC_012024.3	LN:29360087
```


Here we've requested that samtools return only the header section, which contains lots of metadata about the file, including all the contig names in the genome (each @SQ line contains a contig name). Each line begins with an "@" sign. The header can be quite large, especially if there are many contigs in your reference.


## 2.5.  Generating Total Read Counts from Alignment using htseq-count. 
Now we will be using the program htseq-count to count how many reads map to each annotated gene in the genome. To do this, we first need to download the annotation file [here](https://bk-genomica.cebas.csic.es:5001/sharing/A2eYbbxbo). It is in GFF format. More information about the structural annotation (VCost.v3) can be found [here](https://urgi.versailles.inra.fr/Species/Vitis/Annotations)

Once downloaded and unziped, then you can count the features using the htseq-count program.

```
module load htseq/0.11.2
htseq-count  -s no -r pos -f bam ${DIR}/${name}.bam /path/to/Vitis_vinifera_gene_annotation_on_V2_20_myversion.gff3  > ${DIR}/${name}.counts

```

`-s no` indicates we're using an unstranded RNA-seq library.
`-r pos` tells htseq-count that our BAM file is coordinate sorted.
`-f bam` indicates that our input file is in BAM format.
The above command should be repeated for all other BAM files as well. The full script for slurm scheduler can be found in the **count/** folder which is called [htseq_count.sh*().

Once all the bam files have been counted, you can access to all counts [here](https://bk-genomica.cebas.csic.es:5001/sharing/SbV76npgr)
Let's have a look at the contents of a counts file:

head Negroamaro_Touch_1.counts

which will look like:
```
Vitvi06g01793	0
Vitvi06g01792	0
Vitvi06g01791	0
Vitvi06g01790	0
Vitvi06g01796	1
Vitvi06g01795	0
Vitvi06g01794	0
Vitvi13g01117	0
Vitvi13g01116	0
Vitvi13g01115	0
```

We see the layout is quite straightforward, with two columns separated by a tab. The first column gives the vitis gene ID, the second column is the number of mRNA fragments that mapped to the gene. These counts are the raw material for the differential expression analysis in the next section.

## 2.6. Pairwise Differential Expression with Counts in R using DESeq2.

After download the count files to your local computer, to identify differentially expressed (DE) genes, we will use the R package DESeq2, a part of the [Bioconductor](https://www.bioconductor.org/about/) project. After the counts have been generated, typical differential expression analyses can be done easily on laptop computers, instead of on the cluster.
Typical DE analyses involve both a statistical analysis, which is used to rigorously identify genes of interest and their effect sizes, and data visualization, which is used both as a way of illustrating results and as a form of quality control. Even a quick examination of some of the plots we'll make below can reveal problematic samples or unexpected heterogeneity within treatment groups that could bias results, or to put things in a more positive light, inspire new analyses or interpretations.

### Launching `R`
For this tutorial, we'll assume you're using `RStudio`, but however you launch R, it's always a good idea to ensure you begin a new project with a clean workspace. In `RStudio` we recommend you select "File > New Project" and start a new R project in whatever directory you choose.

The following steps will use several different R packages. You'll need to make sure they're installed first.

For differential expression and visualization:

`DESeq2` needs to be installed through Bioconductor. See instructions here.
`apeglm` needs to be installed through Bioconductor. Instructions here
`pheatmap` can be installed by typing `install.packages("pheatmap")`
`ggplot2` is best installed along with the entire set of tidyverse packages. `install.packages("tidyverse")`


### Beginning the statistical analysis
All the following R code has been condensed into a single script [DESeq.R](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/DESeq2/DESeq.R) in the r_analysis directory of the tutorial repository. You can open that in RStudio rather than copying and pasting from here if you'd like.

This first chunk of code will load the necessary R packages and assign values to some R objects we'll use to tell `DESeq2` how and where to read and write files. Note that in this example, the code is expecting to find the "count" directory in the parent of the project directory (that's what directory <- "../count" means). For this code to work, you'll have to edit the path to point to the directory containing your count data.

```
# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")


# create an object with the directory containing your counts:
	# !!edit this to point to your own count file directory!!
directory <- "../count"

# ensure the count files are where you think they are
list.files(directory)

sampleFiles <- list.files(directory, pattern = ".*counts")
```


Now we'll create a data frame that connects sample names, treatments, and count files. `DESeq2` has a function that can use that table to read and format the data so that it can be analyzed.

```
# create a vector of sample names. ensure these are in the same order as the "sampleFiles" object!

sampleNames_white <- c(
  "Garganega_Soft_1.counts_contados",
  "Garganega_Soft_2.counts_contados",
  "Garganega_Soft_3.counts_contados",
  "Garganega_Touch_1.counts_contados",
  "Garganega_Touch_2.counts_contados",
  "Garganega_Touch_3.counts_contados",
  "Glera_Soft_1.counts_contados",
  "Glera_Soft_2.counts_contados",
  "Glera_Soft_3.counts_contados",
  "Glera_Touch_1.counts_contados",
  "Glera_Touch_2.counts_contados",
  "Glera_Touch_3.counts_contados",
  "Moscatobianco_Soft_1.counts_contados",
  "Moscatobianco_Soft_2.counts_contados",
  "Moscatobianco_Soft_3.counts_contados",
  "Moscatobianco_Touch_1.counts_contados",
  "Moscatobianco_Touch_2.counts_contados",
  "Moscatobianco_Touch_3.counts_contados",
  "Passerina_Soft_1.counts_contados",
  "Passerina_Soft_2.counts_contados",
  "Passerina_Soft_3.counts_contados",
  "Passerina_Touch_1.counts_contados",
  "Passerina_Touch_2.counts_contados",
  "Passerina_Touch_3.counts_contados",
  "Vermentino_Soft_1.counts_contados",
  "Vermentino_Soft_2.counts_contados",
  "Vermentino_Soft_3.counts_contados",
  "Vermentino_Touch_1.counts_contados",
  "Vermentino_Touch_2.counts_contados",
  "Vermentino_Touch_3.counts_contados")

# create a vector of conditions. again, mind that they are ordered correctly!
sampleCondition_white <- c("EV",
                     "EV",
                     "EV",
                     "PV",
                     "PV",
                     "PV",
                     "EV",
                     "EV",
                     "EV",
                     "PV",
                     "PV",
                     "PV",
                     "EV",
                     "EV",
                     "EV",
                     "PV",
                     "PV",
                     "PV",
                     "EV",
                     "EV",
                     "EV",
                     "PV",
                     "PV",
                     "PV",
                     "EV",
                     "EV",
                     "EV",
                     "PV",
                     "PV",
                     "PV")

# now create a data frame from these three vectors. 
sampleTable_white <- data.frame(
  sampleName = sampleNames_white,
  fileName = sampleFiles_white,
  condition = sampleCondition_white
)

# look at the data frame to ensure it is what you expect:

head(sampleTable_white)

                         sampleName                          fileName condition
1  Garganega_Soft_1.counts_contados  Garganega_Soft_1.counts_contados        EV
2  Garganega_Soft_2.counts_contados  Garganega_Soft_2.counts_contados        EV
3  Garganega_Soft_3.counts_contados  Garganega_Soft_3.counts_contados        EV
4 Garganega_Touch_1.counts_contados Garganega_Touch_1.counts_contados        PV
5 Garganega_Touch_2.counts_contados Garganega_Touch_2.counts_contados        PV
6 Garganega_Touch_3.counts_contados Garganega_Touch_3.counts_contados        PV


# create the DESeq data object
uva_white <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_white, 
  directory = directory, 
  design = ~ condition 
)
```

Now we have a "DESeqDataSet" object (you can see this by typing `is(ddsHTSeq))`. By default, the creation of this object will order your treatments alphabetically, and the fold changes will be calculated relative to the first factor. If you would like to choose the first factor (e.g. the PV), it's best to set this explicitly in your code. We'll demonstrate that here even though our factor names ("EV" and "PV") already result in the correct order for this data set.

```
# To see the levels as they are now:

ddsHTSeq$condition


 [1] EV EV EV PV PV PV EV EV EV PV PV PV EV EV EV PV PV PV EV EV EV PV PV PV EV EV EV PV PV PV
Levels: EV PV

```

As a next step, we can toss out some genes with overall low expression levels. We won't have statistical power to assess differential expression for those genes anyway. This isn't strictly necessary because `DESeq2` does *independent filtering* of results, but if you have a lot of samples or a complex experimental design, it can speed up the analysis.

```
# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts_uva_white <- rowSums(counts(uva_white_test))
# take the log
logsumcounts_uva_white <- log(sumcounts_uva_white,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts_uva_white,breaks=100)

```
![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/hist_white.png)

```
# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 20
keep_uva_white  <- sumcounts_uva_white  > 20

# keep only the genes for which the vector "keep" is TRUE
ddsHTSeq_uva_white_20filter <- uva_white_test[keep_uva_white,]

```
Ok, now we can do the statistical analysis. All the steps are wrapped in a single function, DESeq().

```
dds_uva_white_filter20 <- DESeq(ddsHTSeq_uva_white_20filter)
```

When we run this function, a few messages will be printed to the screen:
```

estimating size factors
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing

```

You can learn more in depth about the statistical model by starting with [this vignette]() and reading through some of the references to specific parts of the methods, but it will be helpful to explain a little here.



1. *estimating size factors*: This part of the analysis accounts for the fact that standard RNA-seq measures **relative abundances** of transcripts within each sample, not **absolute abundances** (i.e. transcripts/cell). Because of this, if libraries are sequenced to different depths, genes with identical expression levels will have different counts. It may seem that simply adjusting for sequencing depth could account for this issue, but changes in gene expression among samples complicate things, so a slightly more complex normalization is applied.
2. three dispersion estimate steps_: `DESeq2` models the variance in the counts for each sample using a negative binomial distribution. In this distribution the variance in counts is determined by a dispersion parameter. This parameter needs to be estimated for each gene, but for sample sizes typical to RNA-Seq (3-5 per treatment) there isn't a lot of power to estimate it accurately. The key innovation of packages like DESeq2 is to assume genes with similar expression have similar dispersions, and to pool information across them to estimate the dispersion for each gene. These steps accomplish that using an "empirical Bayesian" procedure. Without getting into detail, what this ultimately means is that the dispersion parameter for a particular gene is a compromise between the signal in the data for that gene, and the signal in other, similar genes.
3. fitting model and testing: Once the dispersions are estimated, this part of the analysis asks if there are treatment effects for each gene. In the default case, `DESeq2` uses a *Wald test* to ask if the log2 fold change between two treatments is significantly different than zero.


Once we've done this analysis, we can extract a nice clean table of the results from the messy R object:
```
# get results table
res_uva_white_filter20 <- results(dds_uva_white_filter20)

# get a quick summary of the table
summary(res_uva_white_filter20)

out of 19160 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 6019, 31%
LFC < 0 (down)     : 5791, 30%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


# check out the first few lines
head(res_uva_white_filter20)

log2 fold change (MLE): condition PV vs EV 
Wald test p-value: condition PV vs EV 
DataFrame with 6 rows and 6 columns
                baseMean log2FoldChange     lfcSE      stat      pvalue       padj
               <numeric>      <numeric> <numeric> <numeric>   <numeric>  <numeric>
Vitvi13g01113   2.695885      -0.694770  0.549480 -1.264416 0.206080927 0.28243995
Vitvi14g00358   2.695679       0.647682  0.432861  1.496281 0.134580447 0.19655167
Vitvi14g01916 102.535670       0.183218  0.172438  1.062512 0.288003199 0.37423813
Vitvi14g01910 825.621793       0.661020  0.196417  3.365400 0.000764327 0.00186175
Vitvi14g01911   5.194121       0.114220  0.368994  0.309545 0.756906884 0.81542513
Vitvi14g01918   0.756719       1.794347  0.888438  2.019664 0.043418276 0.07315285

```

In the results table there are six columns:


- basemean: the average of the normalized counts across samples.
- log2FoldChange: the log2-scaled fold change.
- lfcSE: standard error of log2 fold change.
- stat: Wald statistic
- pvalue: raw p-value
- padj: A p-value adjusted for false discoveries.


If we were interested in ranking the genes to look at, say, the top 20 most important, what would be the best factor? Log2 fold changes? Adjusted p-values? This is a tricky question. It turns out that both of these are not great because they are confounded with expression level. Genes with low expression can have extremely inflated fold changes, while genes with very high expression, but low fold changes, can have extremely low p-values.

To deal with this, `DESeq2` provides a function for *shrinking* the log2 fold changes. These shrunken log2 fold changes can be used to more robustly rank genes. Like the empirical Bayesian estimate of the dispersion parameter, this turns the log2 fold changes into a compromise between the signal in the data and a prior distribution that pulls them toward zero. If the signal in the data is weak, the log2 fold change is pulled toward zero, and the gene will be ranked lower. If the signal is strong, it resists the pull. This does not impact the p-values.

We can get the shrunken log fold changes like this:

```

# get shrunken log fold changes
res_shrink_res_uva_white_filter20  <- lfcShrink(dds_uva_white_filter20 ,coef="condition_PV_vs_EV")

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=res_shrink_res_uva_white_filter20$log2FoldChange,
  y=res_shrink_res_uva_white_filter20$log2FoldChange,pch=20,
  cex=.2,
  col=1+(res_shrink_res_uva_white_filter20$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change"
)
abline(0,1)

````

![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/abline.png)


```
# get the top 20 genes by shrunken log2 fold change
uva_white_top20 <- order(-abs(res_shrink_res_uva_white_filter20$log2FoldChange))[1:20]
res_shrink_res_uva_white_filter20[uva_white_top20,]

log2 fold change (MAP): condition PV vs EV 
Wald test p-value: condition PV vs EV 
DataFrame with 20 rows and 5 columns
                 baseMean log2FoldChange     lfcSE       pvalue         padj
                <numeric>      <numeric> <numeric>    <numeric>    <numeric>
Vitvi14g00146     30.1700        11.4745  3.087699  3.50422e-17  3.33868e-16
Vitvi19g00435     48.6401        11.3420  1.112611  9.79171e-09  4.43520e-08
Vitvi15g00947     23.3747        11.1669  2.940035  2.40364e-35  7.12908e-34
Vitvi12g00355   2739.8475       -10.1421  0.801055  5.33556e-38  1.81258e-36
Vitvi05g01947  31191.1907       -10.0558  0.460493 3.93815e-107 2.03932e-104
...                   ...            ...       ...          ...          ...
Vitvi15g00643 10356.93497       -8.77141  0.369831 8.67733e-126 7.55717e-123
Vitvi15g00906     5.76363        8.69918  2.824372  1.16811e-09  5.82533e-09
Vitvi05g00715 32123.94430       -8.67703  0.380777 3.60949e-117 2.46992e-114
Vitvi11g01339     5.90574        8.60915  2.937362  6.44727e-07  2.40799e-06
Vitvi01g02042     4.67364        8.47392  2.679176  1.83052e-16  1.64893e-15



```

### RNA-Seq data visualization
There are several different visualizations we can use to illustrate our results and check to ensure they are robust.

We'll start with the standard **Bland-Altman**, or **MA plot**, which is a high-level summary of the results.

```

plotMA(res_shrink_res_uva_white_filter20, ylim=c(-4,4))

```
![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/plotdispersion.png)

MA-plots depict the log2 fold change in expression against the mean expression for each gene. In this case we're using the shrunken fold changes, but you can plot the raw ones as well. You don't typically learn much from an MA plot, but it does nicely illustrate that you can achieve significance (red dots) at a much smaller fold change for highly expressed genes. It can also be a first indication that a treatment causes primarily up or down-regulation.
    
Next we'll plot the the results of a **principal components analysis (PCA)** of the expression data. For each sample in an RNA-Seq analysis, there are potentially tens of thousands of genes with expression measurements. A PCA attempts to explain the variation in all those genes with a smaller number of principal components. In an RNA-Seq analysis, we typically want to see if the PCs reveal any unexpected similarities or differences among samples in their overall expression profiles. A PCA can often quickly reveal outliers, sample labeling problems, or potentially problematic population structure.

```
# normalized, variance-stabilized transformed counts for visualization
vsd_dds_uva_white_filter20 <- vst(dds_uva_white_filter20, blind=FALSE)

plotPCA(vsd_dds_uva_white_filter20, intgroup="condition")


``` 
![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/PCA.png)


Here we used a *variance stabilized transformation (VST)* of the count data in the PCA, rather than the normalized counts. This is an attempt to deal with the very high range in expression levels (6 orders of magnitude) and the many zeroes, and is similar to, but a bit more robust than simpler approaches, such as simply adding 1 and log2 scaling the normalized counts.

You can see that in our case, the first PC, which explains 70% of the variance in gene expression, separates our two growth stages. In general, most experiments are messier, with much less variance explained. In studies with population structure, the first few PCs often reflect that structure, rather than treatment effects.

*Finally, you have got DE genes between EV and PV in white cultivars (similarly you can get DE genes in red cultivars and in the general analyisis in this study). It is important, due to the that unexpected heterogeneity within replicates of the same genotypes may indicate problems with experimental procedures, or possibly something biologically interesting. Often this will also show up in a PCA or heatmap.*


```
uval2fc_ord <- order(-abs(res_shrink_res_uva_white_filter20$log2FoldChange))
plotCounts(dds_uva_white_filter20, gene=uval2fc_ord[1], intgroup = "condition")

```

![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/gene_lfchihg.png)

Here we have plotted the gene with the largest shrunken log2 fold change in white cultivars.

At the end of the tutorial we'll make a heatmap of the DE genes in white associated with at least one eQTL, 76 in this case (105 in red cultivars and 58 in the general analysis without consider genotypes color). In the plot, genes are clustered by their expression values, and the `vst` expression values are used to color the cells.



**Finally, the obtained counts for down and up regulated genes in red, white and in the general comparsion will be used as input for the downstream eQTL analysis.**

# Variant Calling
In this part we will obtain the genetic variation in the target cultivars. This part is also an introduction to the basics of variant calling from high-throughput, short-read sequencing data. A useful (if dated) review of the underlying concepts is [Nielsen et al. 2011](https://www.nature.com/articles/nrg2986) in Nature Reviews Genetics.

## 3.1 Prepare a reference genome
The first step is to prepare the reference genome. Most software packages that align short-read sequencing data to, or otherwise manipulate a reference genome require that genome to be indexed in some way. We will generate indexes using both `bwa` and `samtools`.


```
#load software modules

module load bwa


module load bwa
module load samtools/1.10

cd ../genome #move to the dir where the genome is hosted.

# set a variable 'GEN' that gives the location of the reference genome:
GEN=pinotnoir.fa.gz

# index the reference genome:
bwa index -p pinotnoir $GEN 

The "-p" is the prefix of the output database.

```

## 3.2 Download data
An overview of the project for the DNA data can be viewed [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA373967) and [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA385116). Genomic DNA from each cultivar was sheared by sonication and 2 × 100-bp paired-end libraries.
The data for each will be two fastq files, one containing the forward reads and one containing reverse reads with members of each pair on the same lines of their corresponding files.

We will use the same 

```
!/bin/bash
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


mkdir ../rawdata #to generate the dir to host the samples
cd ../rawdata #move to this dir

###to start download

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156354/SRR6156354_1.fastq.gz -o SRR6156354_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/004/SRR6156354/SRR6156354_2.fastq.gz -o SRR6156354_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/006/SRR6156356/SRR6156356_1.fastq.gz -o SRR6156356_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR615/006/SRR6156356/SRR6156356_2.fastq.gz -o SRR6156356_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/004/SRR5506714/SRR5506714_1.fastq.gz -o SangioveseC_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/004/SRR5506714/SRR5506714_2.fastq.gz -o SangioveseC_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/005/SRR5506715/SRR5506715_1.fastq.gz -o SangioveseB_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/005/SRR5506715/SRR5506715_2.fastq.gz -o SangioveseB_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/001/SRR5506711/SRR5506711_1.fastq.gz -o SangioveseF_DNA-Seq_of_Vitis_vinifera_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR550/001/SRR5506711/SRR5506711_2.fastq.gz -o SangioveseF_DNA-Seq_of_Vitis_vinifera_2.fastq.gz
.
.
.


```


### Inspecting fastq files

As was commented before in the RNA seq part in a fastq file each read is represented by 4 lines. The first line is the sequence name. The second line is the DNA sequence. The third line, which always begins with a "+" can contain other optional information, but usually does not. The fourth line encodes the base quality scores in ASCII characters. Please go up and check it again if you need a reminder.
Also, it's worth noting that while we are using uncompressed fastq files in this example, all of the programs we are working with will accept files compressed using gzip (generally denoted with the file extension '.gz'), and it's good practice to keep fastq files compressed so they use less storage space.


## 3.3 Asses read quality
In order to evaluate the general quality of reads in the file we will use again the `fastqc` package. Execute this script from the scripts directory by entering `sbatch qualityCheck_Trim.sh` on the command line. The script will do the quality check before and after triming it is a modification from the script used in the RNA seq analysis

Our script runs fastqc as follows:

```

module load sickle
module load fastqc

mkdir ../{fastqc,trimmed} #generate both dir before and after trimming


#Quality check of Reads

DIR=path/to/rawdata

for file in ${DIR}/*_1.fastq.gz
do name=$(basename $file _1.fastq.gz)

echo "========================================== fastqc starting at `date` for  ==========================" $name

fastqc -t 4 -o ../fastqc ${DIR}/${name}_1.fastq.gz  ${DIR}/${name}_2.fastq.gz



```
As previously, once the files are generated you'll have to transfer them to your local computer to open them and examine the results. You can access to the data [here]()


## 3.4 Quality trim
To increase our knowledge, quality trimming is a step in which low quality bases and/or adapter contamination is removed from reads. Current variant callers account for uncertainties in mapping (conditional on the quality of the reference genome) and in base calling, so quality trimming is not always necessary for this application (the worrisome sources of error in variant calling are ["unknown unknowns"](https://en.wikipedia.org/wiki/There_are_known_knowns), like the incompleteness of the reference genome, or systematic error arising from library prep). However, if you have a dataset plagued by adapter contamination or poor quality reads, you may want to try trimming to salvage it and/or remove some of the noise.

For DNA, `sickle` will be the tool used for trimming.

We could trim the reads for the samples as follows:
```
sickle pe \
-t sanger \
-f ${DIR}/${name}_1.fastq.gz \
-r ${DIR}/${name}_2.fastq.gz \
-o ../trimmed/${name}_trimmed_1.fastq \
-p ../trimmed/${name}_trimmed_2.fastq \
-l 45 -q 25 \
-s ../trimmed/singles_${name}.fastq
```

This would discard any read trimmed shorter than 45bp, and if its pair was longer than 45bp, it would be placed in the file given by `-s .../trimmed/singles_${name}.fastq`.

As commented before, our script, which also runs fastqc on the trimmed data, enter `sbatch qualityCheck_Trim.sh` on the command line to run.

## 3.5 Align and compress
Now that we have QC-ed our sequence data, it's time to align it to a reference genome. For that we'll use `bwa`, one of the most widely used short-read aligners. `bwa` implements several alignment methods, but `mem` is best for our application. We previously indexed our reference genome, so we're ready to go here.


```

module load bwa/0.7.17
module load samtools/1.7


cd /open/trimmingdir

GEN=/path/to/genome/pinotnoir



for file in *_1.fastq.gz #loop for all triming files
do name=$(basename $file _trimmed_1.fastq.gz) 

echo "===================bwa -mem started at `date` for ============" $name
bwa mem -t 12  $GEN  ${name}_trimmed_1.fastq.gz ${name}_trimmed_2.fastq.gz  -o ../align/${name}.sam #mem command to align

echo "===================bwa -mem finished at `date` for ============" $name



```

Here -t 12 indicates that the program should use twelve processors, and we feed bwa mem both the location of the reference genome and both files of paired end reads.
Finally, we'll compress the resulting alignment file:
```
#SAM to BAM CONVERSION
echo "===================sam started at `date` for ============" $name

samtools view -bhS ../align/${name}.sam > ../align/${name}.bam #sam conversion to bam

echo "===================sam finished at `date` for ============" $name
```

Execute these scripts from the scripts directory by entering  `sbatch alig.sh` on the command line.

## 3.6 Sort reads by genome position
To call variants at a given position in the reference genome, we need to look at all the reads that overlap that position. In order to do this efficiently, we need to sort the reads in the alignment files by their positions in the reference genome. We'll use `Picard Tools` for this. For example:

```

module load picard/2.9.2

cd ../align #move to the align folder or other dir that hosts the alignment

for file in *.bam #loop for all bam file
do name=$(basename $file .bam) ## name of the sample

java -jar $PICARD SortSam \ #command to run picard
        INPUT=${name}.bam \
        OUTPUT=../align_stepwise/${name}_sort.bam \
        SORT_ORDER=coordinate \
        CREATE_INDEX=True

done
```

From the scripts directory, enter `sbatch Part1e_sort.sh` on the command line.

## 3.7 Mark duplicates
Duplicate sequences are those which originate from the same molecule after extracting and shearing genomic DNA. There are two types: optical and polymerase chain reaction (PCR) duplicates. Optical duplicates are an error introduced by the sequencer. PCR duplicates are introduced by library prepartion protocols that use PCR. Duplicates cause 2 types of artifacts that mislead variant callers.

- First, errors introduced by the polymerase can be propagated to multiple copies of a given fragment. Because these errors are actually part of a DNA sequence, they are likely to have high base qualities. If many sequences from the fragment containing the error are present, the variant caller can be deceived into identifying it as true biological variation.
- Second, when variant callers call genotypes, they assume that heterozygous sites will have equal representation of both alleles in the sequence pool (as they should for germ-line mutations). Dramatically unbalanced coverage of an allele can be a signal that variation is spurious. Because of its exponential reproduction of fragments, PCR can randomly alter allele balance, or amplify small deviations in the initial sample, causing a variant caller to incorrectly call genotypes as homozygotes, or a truly variable site as invariant.

For these reasons we need to exclude duplicate sequences from variant calling. They can be identified most easily from paired-end data as those sequences for which both reads have identical start sites. This may eliminate some sequences which are in fact derived from unique fragments in the original library, but if fragmentation is actually random, identical fragments should be rare. Once identified, duplicate sequences can be marked and ignored during variant calling (or other types of analyses) downstream.

Here is some example code:

```
module load picard/2.9.2
DIR2=path/to/bamfiles


for file in ${DIR2}/*_sort.bam #loop for each bam file in the dir
do name=$(basename $file _sort.bam) #to obtain the name of each sample

java -jar $PICARD MarkDuplicates \ #basic command to mark duplicates
        INPUT=${DIR2}/${name}_sort.bam \
        OUTPUT=${DIR2}/${name}_mkdup.bam \
        REMOVE_DUPLICATES=Ture \
        METRICS_FILE=${DIR2}/${name}_mkdup_metrics.txt \
        CREATE_INDEX=True
done

```
In this example, duplicate sequences remain in the file, but they are flagged as such.

Execute the script from the XX directory by entering `sbatch markduplicates.sh` on the command line.


## 3.8 Index alignment files
The last step in preparing the reads is to index the bam files. This needs to be done to enable fast access to reads overlapping a given position of the genome. Without the index, if you wanted to access reads at the beginning of chromosome 8, you'd need to read through chromosomes 1-7 until you got there. With many samples or deep coverage, this would be a big problem. The bam index is a map to the bam file that lets you skip around quickly. In this case, picard already generated the index during duplicate marking, but if you used a different approach, you'd need to use samtools to do this. For example:

```
module load samtools

# "*mkdup.bam" will refer to each of the 
for file in path/to/markduplicates/*_mkdup.bam
	do samtools index $file
done

```

Now we have completed the initial QC, alignment and processing steps. At this point, you may have noticed that we have accumulated six copies of our data. Two copies of the fastq files, and four copies of the alignment files. This is a large and space-wasting mess. If we were working with many samples of high coverage human genomes, we would want to go and delete the intermediate alignment files and the trimmed fastqs, keeping only the original fastqs and the analysis-ready bams. Another approach, detailed in Part 3, would pipe many of these steps together and avoid creating some of the intermediate files to begin with.

Execute the indexing script from the scripts directory by entering Part1g_indexbams.sh on the command line.

## 3.9 Exploring SAM files
Now we can explore an alignment file. We can get a lot of basic stats on the SAM file using `samtools stats`:

```
module load bcftools
module load samtools

mkdir bamstats

DIR2=path/to/bam



for file in ${DIR2}/*_mkdup.bam
do name=$(basename $file _mkdup.bam)


GEN=/path/to/genome/PinotNoir.fa 

samtools stats -r $GEN ${DIR2}/${name}_mkdup.bam > ${DIR2}/${name}_samstat.txt


```
If you do less in any txt file obtained and you'll see this file is pretty messy, but we can pull out specific parts of it using grep.

These are some basic stats about the alignment file (in this case for cultivar Refosco):

```
grep ^SN Refosco_samstat.txt | cut -f 2- | head -n12
raw total sequences:	164876640
filtered sequences:	0
sequences:	164876640
is sorted:	1
1st fragments:	82438320
last fragments:	82438320
reads mapped:	152540110
reads mapped and paired:	149794714	# paired-end technology bit set + both mates mapped
reads unmapped:	12336530
reads properly paired:	134970986	# proper-pair bit set
reads paired:	164876640	# paired-end technology bit set
reads duplicated:	5012283	# PCR or optical duplicate bit set
```

This is a histogram of per base coverage (showing only 12 lines):

```
grep ^COV Refosco_samstat.txt | cut -f 2- | head -n12
[1-1]	1	3407825
[2-2]	2	3247513
[3-3]	3	3247714
[4-4]	4	3325293
[5-5]	5	3436060
[6-6]	6	3557262
[7-7]	7	3717062
[8-8]	8	3902527
[9-9]	9	4142693
[10-10]	10	4399872
[11-11]	11	4709602
[12-12]	12	5020517


```
And there is much more information. You can access to the samstat files [here](https://bk-genomica.cebas.csic.es:5001/sharing/54YV8bExO)

## 3.10 Variant Discovery method
After the alignment of short read sequencing data to a reference genome and prepared the alignment files for variant calling, we should evaluate the resulting sequence alignments for evidence of variation using `bcftools`, which applies a statistical model to evaluate the evidence for variation.

Why do we apply a statistical model to discover variants? Why can't we just say "this observed sequence differs from the reference genome at this position, therefore we have found a variant"?

In short, because the sequence data we observe have been passed to us through several messy and error-prone laboratory and computational processes. We can break down this pipeline into several steps:

- *Sampling the genome*: A diploid individual has two copies of each chromosome. If this individual is heterozygous for a given site, the probability that we only observe the reference allele is 0.5^n, where n is the number of (non-duplicate) sequences. We still have a 3.125% chance of failing to observe the alternate allele given 5x coverage of the site.
- *Library preparation*: During library preparation, polymerases can misincorporate bases. These bases will be read by the sequencer and likely assigned high base quality scores. If PCR is used during library preparation, random changes in the allele balance at heterozygous sites can be introduced, or existing random biases from the pool of extracted DNA can be amplified, exacerbating the problem in the previous step.
Sequencing: When we load our library of DNA fragments onto the sequencer to be analyzed, the machine may read bases incorrectly.
- *Reference mapping*: After sequencing our library, we must map the sequences back to the reference genome. Most reference genomes are a) incomplete and b) have many copies of similar or identical sequences. These issues can mean that a DNA fragment mapped back to the reference genome may not actually have been derived from that region. Differences between it and the reference therefore don't represent the kind of polymorphism we're looking for.
- *Sequence alignment*: Even if we have assigned a DNA fragment to the correct location in the genome, we still need to align it properly. If there are indels, this can become very difficult and there may be several equally likely alignments.

So, all of this is to say that when we observe a set of sequences aligned to the reference genome, and there appear to be some bases that differ, we need a statistical model that can account for these complex sources of error to help us decide whether that variation is real and assign genotypes to individuals.


## 3.11 Generate a pileup file
`bcftools` uses a two step procedure to call variants. First it generates a pileup file using `bcftools mpileup`. The pileup file summarizes the per base coverage at each site in the genome. Each row represents a genomic position, and each position that has sequencing coverage is present in the file. The second step, using `bcftools call` actually applies the statistical model to evaluate the evidence for variation represented in that summary and generate output in the `Variant Call Format (VCF)`.

Because pileup files for whole genome sequencing contain a summary for every site in the genome for every sample, they can be very large. We also don't typically need to look at them directly, or do anything with them outside of the variant calling. So if you use this approach on your own data, you will usually simply `pipe` the output of bcftools mpileup directly to `bcftools call`.

For teaching purposes here, we'll keep the steps separate.

Here is a call to `bcftools mpileup`:

```
module load bcftools


INDIR=../align_stepwise #dir where the alignments are located

# make output directory if it doesn't exist. 
OUTDIR=../variants_bcftools
mkdir -p $OUTDIR


# make a list of bam files
ls $INDIR/*_mkdup.bam >$INDIR/list.bam

GEN=/path/to/the/genome/PinotNoir.fa 

bcftools mpileup \
        -f $GEN \
        -b $INDIR/list.bam \
        -q 20 -Q 30  > $OUTDIR/vitis.pileup

	
```

We give bcftools the reference genome with -f, a list of bam files with -b, tell it to exclude bases with quality lower than 30 and reads with mapping quality lower than 20 with -q and -Q.


## 3.12 Call variants
The last step is to evaluate the evidence (summarized in the pileup file) that the sequence variation we observe is true biological variation, and not errors introduced during library preparation, sequencing, mapping and alignment. Here we use bcftools call.

```
bcftools call -m -v -Oz -o $INDIR/vitis.vcf.gz $INDIR/vitis.pileup
```

We use the -o flag to indicate the output file name. The flag -m specifies one of two possible variant calling routines, -v says that only variable sites should be output, and -Oz indicates the output should be compressed in a version of the gzip compression format.

In this case, we've told bcftools to output a compressed file. Other variant callers may not have that option. In those cases, it's a good idea to use compression to save space. `bgzip` is a commonly used utility for compressing tabular data in genomics. It has the advantage of producing files that can be indexed, generally using the companion program `tabix`. `bcftools` used the same algorithm as bgzip to compress our file. This indexing is critical for VCF files containing millions of variants.

**Finally, the obtained vcf will we used as input for the downstream eQTL analysis.**

# eQTL analysis
For the analysis of **eQTLs** will be use the R package `Matrix eQTL`, the complete information of this software can be found [here](http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/) and the paper describe this software is [Shabalin 2012](https://pubmed.ncbi.nlm.nih.gov/22492648/). This software is the oficial tool of the [Genotype-Tissue Expression (GTEx) project](https://gtexportal.org/home/). 

This software can accommodate large expression and genotype datasets.`Matrix eQTL` checks for association between each SNP and each transcript by modeling the effect of genotype as either categorical (`ANOVA model`) or additive linear (`least squares model`). The simple linear regression (used in this study) is one of the most commonly used models for **eQTL analysis**. In addition, `Matrix eQTL` can test for the signiﬁcance of genotype-covariate interaction (not considered in this study). `Matrix eQTL` also supports correlated errors to account for relatedness of the samples and heteroscedastic. `Matrix eQTL` implements a separate test for each gene-SNP pair and corrects for multiple comparisons by calculating false discovery rate (`FDR`). 

Five different input ﬁles (*snps*=snps data; gene=expression mean by sample of the normalized read counts of each DEG data; *cvrt*=covariates; *genepos* = gene location; *snpspos* = SNP location) are required to run `Matrix eQTL`. 

To obtain the different expression data we need to obtain the mean value for each sample (average of the 3 replicates) for that a basic command using awk

```
head mean_total_PV.txt
Id	Barbera	Garganega	Glera	Moscatobianco	Negroamaro	Passerina	Primitivo	Refosco	Sangiovese	Vermentino
Vitvi13g01113	1.49213	5.8028	0.331948	2.22721	1.64488	1.35951	6.78414	0.242561	1.88299	0.640516
Vitvi14g00358	1.67492	3.13542	3.39697	3.04996	1.93371	1.34049	6.2225	1.8085	1.63172	5.46181
Vitvi14g01916	78.669	126.816	136.598	114.643	121.646	104.356	115.347	165.344	179.314	64.0109
Vitvi14g01917	0	2.39883	0.642034	0.797722	0.268726	0.503129	1.47842	0	1.58047	0
Vitvi14g01910	853.224	1816.37	793.646	597.777	1299.89	1037.8	836.554	1217.96	770.137	825.26
Vitvi14g01911	8.98841	9.90573	6.19271	4.88006	13.9135	2.17786	14.1546	6.27935	6.01241	3.85567
Vitvi14g01918	0.45301	3.06014	0.973982	0.793992	0	0.84687	0.931662	0.556233	0	0.325113


```


```
head mean_total_EV.txt
Id	Barbera	Garganega	Glera	Moscatobianco	Negroamaro	Passerina	Primitivo	Refosco	Sangiovese	Vermentino
Vitvi13g01113	2.66732	6.74091	2.16952	2.78907	4.47536	3.75126	4.42606	0	1.85146	1.08383
Vitvi14g00358	1.57085	3.62084	3.09695	1.69132	2.57991	0.358892	3.08438	1.53921	4.14435	1.78452
Vitvi14g01916	114.528	132.655	93.551	76.2368	283.154	124.469	105.808	157.033	173.289	50.8079
Vitvi14g01917	0	0	0.37405	0.427394	0	0	0.723733	0.414316	1.04699	0.307503
Vitvi14g01910	1044.96	982.565	643.919	536.382	425.309	362.087	1037.76	313.96	503.027	657.923
Vitvi14g01911	9.04248	6.76973	8.11213	6.48155	39.3775	2.40565	15.412	10.9655	4.29964	1.08885
Vitvi14g01918	1.16305	0.43424	0	0.443603	0.61424	0	0	0	0	0.700687
Vitvi14g01919	40.1727	39.275	66.3442	37.9537	40.1266	42.804	42.3571	43.0221	63.5099	28.9353
Vitvi14g01867	706.947	589.912	451.828	332.496	471.284	477.148	667.631	417.869	366.749	440.258

```


```
head genloc.txt
gene	chr	s1	s2
Vitvi01g00001	chr1	10715	27598
Vitvi01g00002	chr1	31968	36088
Vitvi01g00003	chr1	36772	36945
Vitvi01g00004	chr1	42488	45996
Vitvi01g01833	chr1	46794	47258
Vitvi01g00005	chr1	48568	57979
Vitvi01g01834	chr1	60839	61441
Vitvi01g00006	chr1	93202	107691
Vitvi01g01835	chr1	109284	110282
Vitvi01g01836	chr1	116646	118695

```

```

head meanexpresion_red_EV.txt
gene	Barbera_EV	Negroamaro_EV	Primitivo_EV	Refosco_EV	Sangiovese_EV
Vitvi13g01113	2.66867	4.46885	4.42344	0	1.84415
Vitvi14g00358	1.57186	2.5752	3.08424	1.53795	4.12994
Vitvi14g01916	114.574	282.703	105.775	156.86	172.693
Vitvi14g01910	1045.3	424.632	1037.4	313.625	501.429
Vitvi14g01911	9.04598	39.3132	15.4136	10.9538	4.29156
Vitvi14g01919	40.1823	40.0652	42.3418	42.976	63.2855
```


```
head meanexpresion_red_PV.txt
gene	Barbera_PV	Negroamaro_PV	Primitivo_PV	Refosco_PV	Sangiovese_PV
Vitvi13g01113	1.49397	1.64738	6.78665	0.241977	1.87905
Vitvi14g00358	1.6746	1.93622	6.23241	1.80492	1.62962
Vitvi14g01916	78.6793	121.841	115.449	165.085	179.039
Vitvi14g01910	853.084	1301.89	837.181	1216.18	768.956
Vitvi14g01911	8.98666	13.9366	14.1683	6.2705	6.00358
Vitvi14g01919	117.402	70.9669	61.2112	117.11	114.001

```

```
head meanexpression_EV_white.ttx

Id	Garganega_EV	Glera_EV	Moscatobianco_EV	Passerina_EV	Vermentino_EV
Vitvi13g01113	6.78079	2.17505	2.81288	3.76762	1.08905
Vitvi14g00358	3.63954	3.10557	1.70573	0.36039	1.7914
Vitvi14g01916	133.408	93.8193	76.8836	125.012	50.9717
Vitvi14g01910	988.277	645.772	540.956	363.687	660.145
Vitvi14g01911	6.81011	8.13521	6.53709	2.41645	1.09248
Vitvi14g01918	0.436147	0	0.447396	0	0.702353


```

```
head meanexpression_PV_white.ttx
Id	Garganega_PV	Glera_PV	Moscatobianco_PV	Passerina_PV	Vermentino_PV
Vitvi13g01113	5.7774	0.330887	2.23113	1.35479	0.63926
Vitvi14g00358	3.12526	3.38807	3.05461	1.33523	5.45098
Vitvi14g01916	126.271	136.268	114.842	104	63.881
Vitvi14g01910	1808.79	791.727	598.854	1034.43	823.587
Vitvi14g01911	9.86327	6.17825	4.89039	2.1699	3.84806

```


```
head white_gt_1_0_2_noNA_nohomo.txt
SNP_id	Garganega	Glera	MoscatoBianco	Passerina	Vermentino
NC_012016.3_99_4	2	1	2	2	2
NC_012016.3_5066_442	2	1	2	2	2
NC_012016.3_6391_454	2	1	1	2	1
NC_012016.3_6423_455	2	1	2	2	1
NC_012016.3_6774_468	2	0	2	2	1
NC_012016.3_7139_482	1	2	2	2	1
NC_012016.3_7154_483	2	0	2	2	1
NC_012016.3_7156_484	1	2	2	1	2
NC_012016.3_7160_485	1	2	2	2	2

```


```
head red_gt_1_0_2_noNA_nohomo.txt
SNP_id	Barberanera	NegroAmaro_merge	Zinfandel	Refosco	Sangiovese_merge
NC_012016.3_1528_82	2	0	2	0	2
NC_012016.3_1813_114	2	1	2	2	2
NC_012016.3_1818_115	2	1	2	2	2
NC_012016.3_1824_117	1	2	2	2	2
NC_012016.3_1827_119	2	1	2	2	2
NC_012016.3_1828_120	2	0	1	1	2
NC_012016.3_1861_123	1	0	1	1	1
NC_012016.3_1863_124	2	1	2	1	2
NC_012016.3_1868_125	1	0	1	1	1

```

To remove lines with NA (missing data) a small python script was generated:

```
!/bin/usr/python


file="name of the file"

with open(file) as f, open("newfile", 'w') as newfile:
        for line in f:
                if line.startswith("SNP"):
                        
                        newfile.write(line)

                else:
                        line=line.rstrip("\n").split("\t")
                        if ('NA' not in line):
                                line='\t'.join(line)+'\n'
                                newfile.write(line)
  

```

Also to remove lines with only one genotype clase (so not segregation) a second python script was generated:


```
#!/bin/usr/python


file="pathtofile"

with open(file) as f, open("newfile", 'w') as newfile:
        for line in f:
                if line.startswith("SNP"):
                        
                        newfile.write(line)

                else:
                        line=line.rstrip("\n").split("\t")
                        ocurrences_0 = line.count("0")
                        if ocurrences_0 >= 1 and line.count("1") >= 1 and line.count("2") >= 1:
                                line='\t'.join(line)+'\n'
                                newfile.write(line)

```
Now, all SNP files are cleaned to be run in Matrix eQTL.

```
head reorder_gt_1_0_2_noNA_nohomo.txt 
ind_id	Barberanera	Garganega	Glera	MoscatoBianco	NegroAmaro_merge	Passerina	Primitivo	Refosco	Sangiovese_merge	Vermentino
NC_012016.3_5066_442	2	2	1	2	1	2	1	1	1	2
NC_012016.3_6391_454	2	2	1	1	2	2	2	2	2	1
NC_012016.3_6423_455	2	2	1	2	2	2	2	2	2	1
NC_012016.3_6774_468	2	2	0	2	2	2	2	2	2	1
NC_012016.3_6783_469	2	2	2	2	2	2	2	2	1	2
NC_012016.3_6807_471	2	2	2	2	2	2	2	1	1	2
NC_012016.3_7139_482	2	1	2	2	2	2	1	2	1	1
NC_012016.3_7154_483	1	2	0	2	1	2	2	2	2	1
NC_012016.3_7156_484	2	1	2	2	2	1	1	2	1	2

```


All these ﬁles need to have a speciﬁc format. The columns of all three ﬁles must have matching order, corresponding in each column to a sample and with one gene/SNP/covariate in each row. 

In the case of the genotype ﬁle, if a linear model is used, as in this study, the values must be numerical in this data set. For that reason, `extract.gt` function from the R package `vcfR v1.12` was used to read and extract genotypes from our [VCF](https://bk-genomica.cebas.csic.es:5001/sharing/5QMwx9dUa) ﬁltered ﬁle in numeric format. In addition,  homozyogus SNPs and SNPs with high number of missing genotypes were removed in each case.

The p-value threshold for *cis*-eQTLs (`pvOutputThreshold.cis`) in this study was 1e-8 and the maximum distance at which gene-SNP pair is considered local (`cisDist`) was 1000. 
In our study, covariates were not considered. The location of each gene was obtained from the annotation file described in the RNA part, you can get the file [here](https://bk-genomica.cebas.csic.es:5001/sharing/A2eYbbxbo). As reminder, you can get more information about the structural annotation (VCost.v3) can be found [here](https://urgi.versailles.inra.fr/Species/Vitis/Annotations)
The location of each SNP was obtained from the [VCF](https://bk-genomica.cebas.csic.es:5001/sharing/5QMwx9dUa) ﬁle obtained in the previous section.

To find *cis*-eQTL associations with down and up regulated genes (in red and white cultivars and in general (without consider colour  of cultivarsphenotype)), Matrix eQTL must be run several times as follow
1.	To detect cis eQTLs in down regulated genes in white cultivars (high expression in EV)
2.	To detect cis eQTLs in up regulated genes in white cultivars (high expression in PV)
3.	To detect cis eQTLs in down regulated genes in red cultivars (high expression in EV) 
4.	To detect cis eQTLs in up regulated genes in red cultivars (high expression in PV) 
5.	A general analysis to detect cis eQTLs in down regulated genes (high expression in EV) 
6.	A general analysis to detect cis eQTLs in up regulated genes (high expression in PV)



```
#######################A general analysis to detect cis eQTLs in up regulated genes (high expression in PV)

basedir="/path/to/your/main/dir/"

#input files for matrix eqtl
useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_General = paste(basedir, "gt_1_0_2_noNA_nohomo.txt", sep="");
expression_PV_General = paste(basedir, "mean_total_PV.txt", sep="");
covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. #Setting the threshold to a high value for a large dataset may cause excessively large output files.
#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().
#errorCovariance = numeric();
#errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file #delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with #row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

## Load SNP data
snps_tchgeneral = SlicedData$new();
snps_tchgeneral$fileDelimiter = "\t";      # the TAB character
snps_tchgeneral$fileOmitCharacters = "NA"; # denote missing values;
snps_tchgeneral$fileSkipRows = 1;          # one row of column labels
snps_tchgeneral$fileSkipColumns = 1;       # one column of row labels
snps_tchgeneral$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps_tchgeneral$LoadFile(SNP_file_name_General);


## Load gene expression data
gene_tchgeneral = SlicedData$new();
gene_tchgeneral$fileDelimiter = '\t'; # the TAB character
gene_tchgeneral$fileOmitCharacters = 'NA'; # denote missing values;
gene_tchgeneral$fileSkipRows = 1; # one row of column labels
gene_tchgeneral$fileSkipColumns = 1; # one column of row labels
gene_tchgeneral$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_tchgeneral$LoadFile(expression_PV_General);

#output_file_name = 'eQTL_results_R.txt';



### Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = '\t'; # the TAB character
cvrt$fileOmitCharacters = 'NA'; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
cvrt$fileSliceSize = snps$nCols()+1; # read file in one piece
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");


#Output file name
output_file_name_cis = 'eQTL_cis_results_GENERAL_PV_ULTIMO_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234) #for reproducible results we have set up the seed in R
GENERAL_PV = Matrix_eQTL_main(
  snps = snps_tchgeneral,
  gene = gene_tchgeneral,
  cvrt = cvrt,
  useModel = useModel,
  pvOutputThreshold = 1e-9,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);
  
  
  
 ```
 
 The command start to run and you can observe the first lines here:
 
 ```
Matching data files and location files
20622 of 20622 genes matched
11895933 of 11895933 SNPs matched

Task finished in 7.07 seconds
Reordering genes
Task finished in 67.444 seconds
Processing covariates
Task finished in 0.005 seconds
Processing gene expression data (imputation, residualization)
Task finished in 0.023 seconds
Creating output file(s)
Task finished in 0.079 seconds
Performing eQTL analysis
 0.02% done, 6 cis-eQTLs, 13,507 trans-eQTLs
 0.05% done, 23,375 trans-eQTLs
 0.08% done, 24,397 trans-eQTLs
 0.11% done, 6 cis-eQTLs, 33,680 trans-eQTLs
 
 
 
 ```
 
Until 100% is done it will be searching for *cis*-eQTL and also *trans*-eQTLs (but they are not considered here)
 
 ```
 #now we can get the top eqtls

gene_values= read.table("/Users/pedromartinez/Desktop/TFM_data/analisis_general/mean_total_touch.txt", row.names = 1, header = T)
snps_values= read.table("/Users/pedromartinez/Desktop/TFM_data/analisis_general/gt_1_0_2_noNA_nohomo.txt", row.names = 1, header = T)

cis_eqtl_res_GENERAL_touch = GENERAL_touch$cis$eqtls
cis_eqtl_res_GENERAL_touch = cis_eqtl_res_GENERAL_touch[cis_eqtl_res_GENERAL_touch$FDR < 0.1,]
top_eqtls_GENERAL_touch = cis_eqtl_res_GENERAL_touch[order(cis_eqtl_res_GENERAL_touch$pvalue),]
top_eqtls_GENERAL_touch = top_eqtls_GENERAL_touch[!duplicated(top_eqtls_GENERAL_touch$gene),]
mafs = apply(as.matrix(snps_values),1,mean)/2
mafs = data.frame(snps=names(mafs), maf = mafs)
top_eqtls_GENERAL_touch = merge(top_eqtls_GENERAL_touch, mafs, by="snps")
top_eqtls_GENERAL_touch = top_eqtls_GENERAL_touch[order(top_eqtls_GENERAL_touch$FDR),]
head(top_eqtls_GENERAL_touch)

#we can plot a gene of interest and specific SNP
gene_id_GENERAL_PV = "Vitvi08g00849"
snp_id_GENERAL_PV = "NC_012014.3_10631405_15825489"



# generate the dataset
data_GENERAL_PV = data.frame(t(snps_values[snp_id_GENERAL_touch,]), t(gene_values[gene_id_GENERAL_touch,]))
# Get reference and alternative allele of the SNP
ref_alt_GENERAL_PV = unlist(snppos[snppos$CHR_POS_0 == snp_id_GENERAL_touch, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_GENERAL_PV= c(paste(ref_alt_GENERAL_touch[1], ref_alt_GENERAL_touch[1], sep="/"), paste(ref_alt_GENERAL_touch[1],
                                                           ref_alt_GENERAL_touch[2], sep="/"), paste(ref_alt_GENERAL_touch[2], ref_alt_GENERAL_touch[2], sep="/"))
gt_states_GENERAL_touch = factor(gt_states_GENERAL_touch, levels=gt_states_GENERAL_touch)
# Assign the labels
data_GENERAL_touch$gt = gt_states_GENERAL_touch[data_GENERAL_touch[,snp_id_GENERAL_touch]+1]
# Subset to only genotype labels and expression
data_GENERAL_touch = data_GENERAL_touch[,c("gt", gene_id_GENERAL_touch)]
colnames(data_GENERAL_touch) = c("genotype", "expression")

```
The small data set for with the expression of selected gene (Vitvi08g00849) and the genotypes of each sample is:

```
data_GENERAL_PV 
                 genotype expression
Barberanera           A/A   0.000000
Garganega             A/A   0.000000
Glera                 A/A   0.000000
MoscatoBianco         A/A   0.000000
NegroAmaro_merge      G/G   4.516550
Passerina             A/A   0.000000
Primitivo             A/A   0.000000
Refosco               A/A   0.000000
Sangiovese_merge      A/G   4.187960
Vermentino            A/A   0.393185
```
We can use ggplot2 to plot the cis-eQTL.

```

# Plot
library(ggplot2)

plotgeneral_PV <- ggplot(data_GENERAL_touch, aes(genotype, expression))+
 geom_point(position = jitter, alpha=0.9, fill="black", size=5) + theme() +theme(text = element_text(size = 16, color ="black"), axis.text= element_text( colour="black", size=14))

print(plotgeneral_PV)

```
![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/genotype3classes.png)

The complete R script to obtain the results of this study can be found [here](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/DESeq2/DESeq.R). After running the script will have the information for all down and up regulated genes in red and white cultivars and in the general analysis as the previous example.

Finally, a a heatmap of up and dow regulated genes in the general analysis (58 in the general analysis without consider genotypes color) can be observed here. In the plot, genes are clustered by their expression values, and the `vst` expression values are used to color the cells.
The function used was `pheatmap`.

```
genesgeneral <- c("Vitvi12g02290",	
                  "Vitvi12g01431",	
                  "Vitvi14g02682",	
                  "Vitvi14g01945",	
                  "Vitvi15g00038",	
                  "Vitvi15g01253",	
                  "Vitvi18g02506",	
                  "Vitvi18g01186",	
                  "Vitvi03g00486",	
                  "Vitvi11g01420",	
                  "Vitvi13g00818",	
                  "Vitvi01g02015",	
                  "Vitvi06g01130",	
                  "Vitvi06g01133",	
                  "Vitvi14g02924",	
                  "Vitvi13g00819",	
                  "Vitvi19g01851",	
                  "Vitvi15g01495",	
                  "Vitvi14g01131",	
                  "Vitvi04g02151",	
                  "Vitvi10g01636",	
                  "Vitvi06g00071",	
                  "Vitvi09g02012",	
                  "Vitvi18g02648",	
                  "Vitvi02g00189",	
                  "Vitvi01g01902",	
                  "Vitvi13g02109",	
                  "Vitvi19g00645",	
                  "Vitvi14g00205",	
                  "Vitvi16g01497",	
                  "Vitvi18g00384",	
                  "Vitvi03g00703",	
                  "Vitvi10g01290",	
                  "Vitvi11g00497",	
                  "Vitvi07g00006",	
                  "Vitvi16g00122",	
                  "Vitvi09g01367",	
                  "Vitvi03g01373",	
                  "Vitvi15g00736",	
                  "Vitvi16g00471",	
                  "Vitvi18g00648",	
                  "Vitvi07g02192",	
                  "Vitvi01g01552",	
                  "Vitvi02g01005",	
                  "Vitvi04g01804",	
                  "Vitvi07g00301",	
                  "Vitvi02g00317",	
                  "Vitvi19g01851",	
                  "Vitvi08g02058",	
                  "Vitvi14g02924",	
                  "Vitvi02g00189",	
                  "Vitvi14g01131",	
                  "Vitvi14g03045",	
                  "Vitvi18g02879",	
                  "Vitvi18g00917",	
                  "Vitvi12g02656",	
                  "Vitvi02g00196",	
                  "Vitvi13g02109",	
                  "Vitvi19g00645",	
                  "Vitvi05g00510",	
                  "Vitvi18g02715",	
                  "Vitvi05g01202",	
                  "Vitvi08g01611",	
                  "Vitvi01g01992")


general1<- subset(vsd_dds_uva_analisis_general_20filter, rownames(vsd_dds_uva_analisis_general_20filter) %in% genesgeneral)
head(general1)
generalpng<- pheatmap(
  assay(general1), 
  cluster_rows=TRUE, 
  fontsize = 6,
  show_rownames=TRUE,
  cluster_cols=FALSE,
  cellheight = 5,
  annotation_col=df
)

```


![Screenshot](https://github.com/pjmartinez/TFM_UM-eqtls/blob/main/uvageneral.png)


# GO biomart
```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
library(dplyr)

##############################
# Select a mart and dataset
##############################

# see a list of "marts" available at host "plant.ensembl.org"
listMarts(host="plants.ensembl.org")

# create an object for the plant Ensembl Genes
mart <- useMart(biomart="plants_mart", host="plants.ensembl.org")

# occasionally ensembl will have connectivity issues. we can try an alternative function:
# select a mirror: 'www', 'uswest', 'useast', 'asia'
# mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", mirror = "useast")

# see a list of datasets within the mart
# at the time of writing, there were 118
head(listDatasets(mart))
```

               dataset                                 description         version
1      aalpina_eg_gene           Arabis alpina genes (A_alpina_V4)     A_alpina_V4
2   achinensis_eg_gene Actinidia chinensis genes (Red5_PS1_1.69.0) Red5_PS1_1.69.0
3     acomosus_eg_gene                 Ananas comosus genes (F153)            F153
4     ahalleri_eg_gene         Arabidopsis halleri genes (Ahal2.2)         Ahal2.2
5      alyrata_eg_gene            Arabidopsis lyrata genes (v.1.0)           v.1.0
6 aofficinalis_eg_gene      Asparagus officinalis genes (Aspof.V1)        Aspof.V1

```

# figure out which dataset is the croaker
# be careful using grep like this. verify the match is what you want
searchDatasets(mart,pattern="vvinifera_eg_gene")
```

 dataset                description version
117 vvinifera_eg_gene Vitis vinifera genes (12X)     12X

```
# there's only one match, get the name
vitisdata <- searchDatasets(mart,pattern="vvinifera_eg_gene")[,1]


# create an object for the vitis dataset
vitis_mart <- useMart(biomart="plants_mart", host="plants.ensembl.org", dataset = vitisdata)

#########################
# Query the mart/dataset
#########################

# filters, attributes and values

# see a list of all "filters" available for the lcrocea dataset.
# at the time of writing, over 205
head(listFilters(vitis_mart))

```

                name                            description
1    chromosome_name               Chromosome/scaffold name
2              start                                  Start
3                end                                    End
4             strand                                 Strand
5 chromosomal_region e.g. 1:100:10000:-1, 1:100000:200000:1
6          with_embl With European Nucleotide Archive ID(s)

```

# see a list of all "attributes" available
# 122 available at the time of writing
head(listAttributes(mart = vitis_mart, page="feature_page"))
```


                   name              description         page
1       ensembl_gene_id           Gene stable ID feature_page
2 ensembl_transcript_id     Transcript stable ID feature_page
3    ensembl_peptide_id        Protein stable ID feature_page
4       ensembl_exon_id           Exon stable ID feature_page
5           description         Gene description feature_page
6       chromosome_name Chromosome/scaffold name feature_page

```

# we can also search the attributes and filters
searchAttributes(mart = vitis_mart, pattern = "ensembl_gene_id")

```
searchAttributes(mart = vitis_mart, pattern = "ensembl_gene_id")
                name    description         page
1    ensembl_gene_id Gene stable ID feature_page
123  ensembl_gene_id Gene stable ID    structure
157  ensembl_gene_id Gene stable ID     homologs
1308 ensembl_gene_id Gene stable ID          snp
1360 ensembl_gene_id Gene stable ID    sequences
```

searchFilters(mart = vitis_mart, pattern="ensembl")
```

                    name                                       description
33       ensembl_gene_id          Gene stable ID(s) [e.g. ENSRNA049441947]
34 ensembl_transcript_id Transcript stable ID(s) [e.g. ENSRNA049441947-T1]
35    ensembl_peptide_id Protein stable ID(s) [e.g. VIT_00s0120g00010.t01]
36       ensembl_exon_id              Exon ID(s) [e.g. ENSRNA049441947-E1]
```


GD <-c("LOC104877553", #list of genes Down regulated in the general analysis
"VIT_15s0048g00620",
       "VIT_09s0054g01520",
       "VIT_18s0001g08990",
       "VIT_13s0106g00060",
       "VIT_19s0015g00060",
       "VIT_18s0001g05300",
       "VIT_03s0088g00140",
       "VIT_15s0048g00840",
       "VIT_16s0013g01880",
       "LOC100271994",
       "VIT_01s0010g02600",
       "VIT_02s0033g00310",
       "LOC100252434",
       "VIT_07s0005g00170",
       "VIT_02s0025g03450",
       "VIT_18s0086g00500",
       "VIT_18s0001g11930",
       "VIT_05s0020g03530",
       "VIT_01s0011g05620")
 
# get gene names and transcript lengths when they exist
annGD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=GD,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annGD <- group_by(annGD, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

#get GO annotation
go_annGD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=GD,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

#check 
head(go_annGD)
```
    ensembl_gene_id description      go_id                                  name_1006     namespace_1003
1 VIT_01s0010g02600             GO:0005634                                    nucleus cellular_component
2 VIT_01s0010g02600             GO:1990275                        preribosome binding molecular_function
3 VIT_01s0010g02600             GO:0051973 positive regulation of telomerase activity biological_process
4 VIT_01s0010g02600             GO:0042254                        ribosome biogenesis biological_process
5 VIT_01s0010g02600             GO:0016887                                     ATPase molecular_function
6 VIT_01s0010g02600             GO:0005524                                ATP binding molecular_function
```

#write the results
write.table(go_annGD, file="goGD.txt")

dim(go_annGD)

```
[1] 58  5
```

GU <- c("VIT_08s0105g00530",
        "VIT_14s0068g00050",
        "VIT_02s0025g02160",
        "VIT_14s0066g02110",
        "GSVIVT00026525001",
        "VIT_02s0025g02220",
        "VIT_18s0001g11470",
        "VIT_05s0029g01220",
        "VIT_08s0007g04860",
        "VIT_12s0028g02930",
        "VIT_12s0055g01170",
        "VIT_14s0081g00360",
        "VIT_14s0108g00620",
        "VIT_15s0024g01070",
        "VIT_18s0001g01200",
        "VIT_03s0063g00650",
        "LOC104880652",
        "VIT_13s0106g00190",
        "VIT_01s0011g06670",
        "VIT_06s0009g02350",
        "VIT_06s0009g02390",
        "VIT_13s0106g00200",
        "LOC104879152",
        "VIT_10s0116g01580",
        "VIT_06s0004g00690",
        "VIT_01s0011g02420",
        "BQ798972",
        "VIT_16s0039g00930",
        "VIT_11s0016g05640",
        "VIT_16s0039g01910",
        "VIT_03s0038g01520"
)

annGU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=GU,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annGU <- group_by(annGU, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

go_annGU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=GU,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

head(go_annGU)
write.table(go_annGU, file="goGU.txt")
dim(go_annGU)




RU <- c(
        "VIT_11s0016g01000",
        "VIT_11s0016g05640",
        "VIT_11s0016g05650",
        "VIT_11s0052g01320",
        "LOC104880822",
        "VIT_12s0028g02950",
        "VIT_12s0028g03350",
        "VIT_12s0057g00820",
        "TC62812",
        "VIT_13s0139g00420",
        "VIT_14s0060g01650",
        "VIT_14s0060g01830",
        "VIT_14s0030g00230",
        "VIT_14s0030g01990",
        "VIT_14s0081g00360",
        "VIT_14s0006g01780",
        "VIT_14s0083g00250",
        "CF415548",
        "VIT_14s0083g01120",
        "VIT_14s0068g00150",
        "VIT_15s0024g01260",
        "VIT_15s0021g01080",
        "VIT_15s0021g01600",
        "VIT_15s0048g02820",
        "VIT_15s0048g02910",
        "VIT_15s0046g00310",
        "VIT_15s0046g02080",
        "VIT_16s0022g01210",
        "VIT_16s0050g02090",
        "VIT_17s0000g04820",
        "VIT_18s0001g00160",
        "VIT_18s0001g03170",
        "VIT_18s0001g03510",
        
        "VIT_18s0001g06300",
        "VIT_18s0001g07340",
        "VIT_18s0001g08450",
        "VIT_18s0001g11070",
        "VIT_18s0001g11390",
        "VIT_19s0015g01460",
        "VIT_01s0011g06590",
        "VIT_01s0026g00080",
        "VIT_01s0010g00340",
        "LOC104880525",
        "VIT_02s0025g02130",
        "VIT_02s0025g02980",
        "VIT_02s0087g00190",
        "VIT_03s0038g03590",
        "VIT_03s0063g01710",
        
        "VIT_03s0088g00820",
        "VIT_04s0008g07190",
        "VIT_04s0069g00060",
        "VIT_04s0023g01780",
        "LOC104879152",
        "VIT_06s0004g02690",
        "VIT_06s0004g03520",
        "VIT_06s0004g05060",
        "VIT_06s0004g06950",
        "VIT_06s0004g07460",
        "VIT_00s0179g00360",
        "VIT_07s0104g00920",
        "VIT_07s0104g00930",
        "VIT_07s0255g00140",
        "VIT_08s0058g00660",
        "VIT_08s0040g01070",
        "VIT_08s0007g03790",
        "VIT_08s0007g06860",
        "VIT_09s0070g00410",
        "VIT_12s0028g01710",
        "VIT_05s0062g01080",
        "VIT_09s0002g05270")


annRU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=RU,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annRU <- group_by(annRU, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

go_annRU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=RU,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

head(go_annRU)
write.table(go_annRU, file="goRU.txt")
dim(go_annRU)




RD <- c("TC65419",
        "VIT_12s0059g00430",
        "VIT_12s0034g00770",
        "VIT_13s0067g01250",
        "VIT_14s0128g00570",
        "VIT_14s0108g00520",
        "VIT_16s0013g01880",
        "VIT_17s0000g00800",
        "LOC100255295",
        "VIT_18s0001g01480",
        "VIT_18s0001g13670",
        "VIT_18s0086g00500",
        "LOC104882664",
        "VIT_18s0041g00150",
        "LOC100271994",
        "VIT_01s0011g02390",
        "VIT_02s0025g03450",
        "VIT_03s0017g01440",
        "VIT_05s0077g00870",
        "VIT_05s0020g01010",
        "VIT_05s0049g01210",
        "VIT_07s0005g06400",
        "VIT_08s0040g00130",
        "VIT_09s0054g00560",
        "VIT_14s0006g01210",
        "VIT_01s0011g02290",
        "LOC104881810")





annRD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=RD,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annRD <- group_by(annRD, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

go_annRD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=RD,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

head(go_annRD)
write.table(go_annRD, file="goRD.txt")
dim(go_annRD)

WU <- c("VIT_00s0188g00050",
        "VIT_11s0016g01030",
        "VIT_11s0016g05640",
        "VIT_12s0028g01710",
        "LOC100250594",
        "LOC104880822",
        "VIT_12s0034g01340",
        "VIT_13s0067g02030",
        "VIT_14s0060g01800",
        "VIT_14s0083g00250",
        "CF415548",
        "VIT_14s0108g00620",
        "VIT_14s0108g00700",
        "VIT_15s0024g01070",
        "VIT_15s0024g01260",
        "VIT_15s0021g01750",
        "VIT_15s0048g01570",
        "VIT_15s0046g01960",
        "VIT_16s0039g01540",
        "GSVIVT00014257001",
        "VIT_17s0000g05980",
        "VIT_18s0001g01200",
        "VIT_18s0001g03160",
        "VIT_18s0001g15470",
        "VIT_18s0166g00080",
        "VIT_18s0001g00920",
        "VIT_19s0014g01970",
        "VIT_19s0014g02000",
        "VIT_19s0090g01050",
        "VIT_19s0015g02070",
        "VIT_01s0011g00190",
        "GSVIVT00030120001",
        "VIT_01s0011g06670",
        "LOC104879548",
        "VIT_02s0025g01290",
        "VIT_03s0038g03080",
        "VIT_03s0063g00650",
        "VIT_03s0088g00820",
        "VIT_03s0088g01240",
        "VIT_04s0008g01810",
        "VIT_04s0008g01830",
        "VIT_04s0008g07170",
        "VIT_04s0069g00060",
        "LOC104879152",
        "VIT_05s0020g00980",
        "VIT_05s0020g01030",
        "VIT_05s0020g01350",
        "VIT_05s0020g04220",
        "VIT_06s0004g03520",
        "VIT_06s0009g00410",
        "VIT_06s0009g02640",
        "VIT_07s0005g01460",
        "VIT_00s0779g00020",
        "VIT_08s0032g00510",
        "VIT_08s0040g01660",
        "VIT_08s0007g02450",
        "VIT_08s0007g03790",
        "VIT_09s0002g01350",
        "VIT_18s0001g06300",
        "VIT_04s0008g01090")

annWU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=WU,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annWU <- group_by(annWU, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

go_annWU <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=WU,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

head(go_annWU)
write.table(go_annWU, file="goWU.txt")
dim(go_annWU)



WD <- c("VIT_14s0060g02020",
        "VIT_14s0128g00340",
        "VIT_17s0000g09720",
        "VIT_18s0001g09660",
        "LOC104882791",
        "LOC100266415",
        "VIT_03s0038g03450",
        "VIT_08s0105g00370",
        "VIT_09s0002g02800",
        "VIT_09s0054g00560",
        "VIT_18s0086g00320")



annWD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=WD,attributes=c("ensembl_gene_id","description","transcript_length"),mart=vitis_mart)


# pick only the longest transcript for each gene ID
annWD <- group_by(annWD, ensembl_gene_id) %>% 
  summarize(.,description=unique(description),transcript_length=max(transcript_length))

go_annWD <- getBM(useCache = FALSE,filter="ensembl_gene_id",value=WD,attributes=c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),mart=vitis_mart)

head(go_annWD)
write.table(go_annWD, file="goWD.txt")
dim(go_annWD)



```
We have annotated all the DEGs with at least one GO term.

**In summary, this tutorial will help you to integrate gene expression information with genome structural features to identify genes and relevant mechanisms for complex processes. Here we have focused in the understanding of the ripening of fruits of Vitis vinifera L. An inventory of *cis*-eQTLs will be an important resource for future research to understand the mechanism for variation in gene regulation during ripening in this species, and could be considered general markers of ripening in grapes. For other woody species the same approach can be used when the availability of such comprehensive data sets becomes a reality.**
