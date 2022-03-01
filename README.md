# TFM_UM-eqtls
## Title:Bioinformatic approach for the discovery of cis-eQTL signals during fruit ripening of a woody species as the grape (Vitis vinifera L.).


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a [Slurm](https://slurm.schedmd.com/sbatch.html) scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. More information about [Slurm](https://slurm.schedmd.com/sbatch.html) can be found in the [hpc wiki](https://hpc-wiki.info/hpc/SLURM). Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
In summary, the repository includs the different inputs, scripts, outputs files (the supplemental files) obtained during the analysis performance in the biorxv paper: https://www.biorxiv.org/content/10.1101/2021.07.06.450811v1





![Screenshot|width=150px](/Figures/generalpipeline.png)

## Cloning the tutorial repository

To work through this tutorial, copy it to your home directory using the git clone command:
```
git clone < git-repository-path.git > 

```

Contents
 1. [Overview](#-1-Overwiew)
 2. [RNA Seq](#-2-RNA-Seq-Analysis)
 3. [Variant Calling](#-3-Variant-Calling)
 4. [eQTL analysis](#-4-eQTL-analysis)

# 1. Overview
This tutorial will teach you how to use open source quality control, RNA Seq, Variant Calling, eQTL tools to complete a cis-eQTL analysis which is possible when you  have generated the specific datasets. Moving through the tutorial, you will take expression and genotypic data from a woody species as grape and perform a eQTL analysis via Matrix eQTL to characterize the gene expression levels during ripening Vitis vinifera L. fruit.

# 2. RNA Seq Analysis
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
#SBATCH --job-name=fastq_dump_xanadu
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=15G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
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

index/
|-- pinotnoir.1.ht2
|-- pinotnoir.2.ht2
|-- pinotnoir.3.ht2
|-- pinotnoir.4.ht2
|-- pinotnoir.5.ht2
|-- pinotnoir.6.ht2
|-- pinotnoir.7.ht2
`-- pinotnoir.8.ht2

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

module load samtools/1.9
samtools view -H LB2A_SRR1964642.bam | head
which will give:

@HD	VN:1.0	SO:coordinate
@SQ	SN:NC_040011.1	LN:43682218
@SQ	SN:NC_040012.1	LN:14376772
@SQ	SN:NC_040013.1	LN:52095323
@SQ	SN:NC_040014.1	LN:6444570
@SQ	SN:NC_040015.1	LN:5657075
@SQ	SN:NC_040016.1	LN:27037660
@SQ	SN:NC_040017.1	LN:29365971
@SQ	SN:NC_040018.1	LN:33955600
@SQ	SN:NC_040019.1	LN:13800884



Here we've requested that samtools return only the header section, which contains lots of metadata about the file, including all the contig names in the genome (each @SQ line contains a contig name). Each line begins with an "@" sign. The header can be quite large, especially if there are many contigs in your reference.

We can use samtools to access reads mapping to any part of the genome using the view submodule like this:

samtools view LB2A_SRR1964642.bam NC_040019.1:171000-172000

This will print to the screen all the reads that mapped to the genomic interval NC_040019.1:171000-172000.

You can use pipes and other linux tools to get basic information about these reads:

samtools view LB2A_SRR1964642.bam NC_040019.1:171000-172000 | wc -l

wc -l counts lines of text, so this command indicates that 411 reads map to this 1kb interval.




## 2.5.  the function htseq-count from the HTSeq v0.13.5 package was used to count how many reads map to each annotated exon (gene) in the genome. 
Now we will be using the program htseq-count to count how many reads map to each annotated gene in the genome. To do this, we first need to download the annotation file. It is in GFF format. It can be done using the following command:

wget ftp://ftp.ensembl.org/pub/release-104/gtf/larimichthys_crocea/Larimichthys_crocea.L_crocea_2.0.104.gtf.gz
gunzip Larimichthys_crocea.L_crocea_2.0.104.gtf.gz
Once downloaded and unziped, then you can count the features using the htseq-count program.

htseq-count -s no -r pos -f bam ../align/LB2A_SRR1964642.bam Larimichthys_crocea.L_crocea_2.0.104.gtf > LB2A_SRR1964642.counts
-s no indicates we're using an unstranded RNA-seq library.
-r pos tells htseq-count that our BAM file is coordinate sorted.
-f bam indicates that our input file is in BAM format.
The above command should be repeated for all other BAM files as well. The full script for slurm scheduler can be found in the count/ folder which is called htseq_count.sh.

Once all the bam files have been counted, the following files will be found in the count directory.

count/
├── htseq_count_NNNNN.err
├── htseq_count_NNNNN.out
├── htseq_count.sh
├── Larimichthys_crocea.L_crocea_2.0.104.gtf
├── LB2A_SRR1964642.counts
├── LB2A_SRR1964643.counts
├── LC2A_SRR1964644.counts
└── LC2A_SRR1964645.counts
Let's have a look at the contents of a counts file:

head Negroamaro_Touch_1.counts

which will look like:

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


We see the layout is quite straightforward, with two columns separated by a tab. The first column gives the Ensembl gene ID, the second column is the number of mRNA fragments that mapped to the gene. These counts are the raw material for the differential expression analysis in the next section.





## 2.6. These ﬁnal counts per gene are the inputs of the R package DESeq2 v3.13, used for the differential expression analysis.


# 3. Variant Calling

# 4. eQTL analysis
