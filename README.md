# TFM_UM-eqtls
## Title: Algorithms for the discovery of cis-eQTL signals in woody species: the vine (*Vitis vinifera* L.) as a study model.


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a [Slurm](https://slurm.schedmd.com/sbatch.html) scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. More information about [Slurm] can be found in the [hpc wiki](https://hpc-wiki.info/hpc/SLURM). Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
In summary, the repository includs the different inputs, scripts, outputs files (the supplemental files) obtained during the analysis performance in the biorxv paper: https://www.biorxiv.org/content/10.1101/2021.07.06.450811v1





![Screenshot](/Figures/generalpipeline.png)


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

The full script for slurm scheduler for red and white cultivars can be found in the [raw_data](TFM_UM-eqtls/rawdata/) folder. Before running it, add your own e-mail address to the --mail-user option (or delete the line entirely if you don't want an e-mail notification when the job completes).

When you're ready, you can execute the script by entering sbatch [rna_red_download.sh](TFM_UM-eqtls/rawdata/rna_red_download.sh) or [rna_white_download.sh](TFM_UM-eqtls/rawdata/) in the terminal. This submits the job to the SLURM scheduler.

Once the job is completed the folder structure will look like this :

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
In total 120 files should be there. The sequence files are in fastq format and compressed using gzip (indicated by the .gz). It's good practice to keep sequence files compressed. Most bioinformatics programs can read them directly without needing to decompress them first, and it doesn't get in the way of inspecting them either. The .out and .err files are output produced by SLURM that you can use to troubleshoot if things go wrong. Lets have a look at at the contents of one of the fastq files:

zcat SRR1631842_GSM1532795_Barbera_Harv_3_Vitis_vinifera_RNA-Seq.fastq.gz | head -n 12
```
@SRR1631844.1 FCC355RACXX:2:1101:1476:2162 length=90
CAACATCTCAGTAGAAGGCGGCGCCTTCACCTTCGACGTGGGGAATCGCTTCAACCTCACGGGGGCTTTCCTCTACACGTCCTGTCCGGA
+SRR1631844.1 FCC355RACXX:2:1101:1476:2162 length=90
?@@D?DDBFHHFFGIFBBAFG:DGHDFHGHIIIIC=D<:?BBCCCCCBB@BBCCCB?CCBB<@BCCCAACCCCC>>@?@88?BCACCBB>
@SRR1631844.2 FCC355RACXX:2:1101:1641:2127 length=90
NGCCTGTAAAATCAAGGCATCCCCTCTCTTCATGCACCTCCTGAAATAAAAGGGCCTGAATAATGTCGTACAGAAGACTGCGGCACAGAC
+SRR1631844.2 FCC355RACXX:2:1101:1641:2127 length=90
#1=DDFFFHHHHGJJJJJIIIJIJGIIJJJIJIJJGIJIJJJJIJJJJJJIJJJIJJJJJJJGIIHIGGHHHHHFFFFFDEDBDBDDDDD
@SRR1631844.3 FCC355RACXX:2:1101:1505:2188 length=90
GGACAACGCCTGGACTCTGGTTGGTATTGTCTCCTGGGGAAGCAGCCGTTGCTCCACCTCCACTCCTGGTGTCTATGCCCGTGTCACCGA
+SRR1631844.3 FCC355RACXX:2:1101:1505:2188 length=90
CCCFFFFFHHFFHJJJIIIJHHJJHHJJIJIIIJEHJIJDIJJIIJJIGIIIIJGHHHHFFFFFEEEEECDDDDEDEDDDDDDDADDDDD
```

Each sequence record has four lines. The first is the sequence name, beginning with @. The second is the nucleotide sequence. The third is a comment line, beginning with +, and which here only contains the sequence name again (it is often empty). The fourth are phred-scaled base quality scores, encoded by ASCII characters. Follow the links to learn more, but in short, the quality scores give the probability a called base is incorrect.



## 2.2. Quality Control Using Trimmomatic
`Trimmomatic` is commonly used to trim low quality and adapter contaminated sequences.

Our usage looks like this for a single sample:
```

module load Trimmomatic/0.39

DIR=path to raw_data
DIR2= path to adapters

for file in ${DIR}/*.fastq.gz 
do name=$(basename $file _Vitis_vinifera_RNA-Seq.fastq.gz)

java -jar  /home/cebas/pmartinez/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        ${DIR}/${name}_Vitis_vinifera_RNA-Seq.fastq.gz  \
        ${DIR}/${name}_trim.fastq.gz \
        ILLUMINACLIP:${DIR2}/TruSeq2-SE.fa:2:30:10 \
        SLIDINGWINDOW:4:20 \
        MINLEN:45

```
 
We call SE for single-end mode and we specify the input and output file names. The ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 command searches for adapter sequence, so we provide a fasta file containing the adapters used in the library preparation, and the numbers control the parameters of adapter matching (see the manual for more details). SLIDINGWINDOW:4:20 scans through the read, cutting the read when the average base quality in a 4 base window drops below 20. We linked to an explanation of phred-scaled quality scores above, but for reference, scores of 10 and 20 correspond to base call error probabilities of 0.1 and 0.01, respectively. MINLEN:45 causes reads to be dropped if they have been trimmed to less than 45bp.Here is a useful paper on setting trimming parameters for RNA-seq.

The full scripts for slurm scheduler is calling which can be found in the quality_control/ folder. Navigate there and run the script by enteriing sbatch fastq_trimming.sh on the command-line.

Following the trimmomatic run, the resulting file structure will look as follows:

quality_control/





Examine the .out file generated during the run. Summaries of how many reads were retained for each file were written there. Here's one example:

TrimmomaticSE: Started with arguments:
 -threads 12 ../raw_data/LB2A_SRR1964642.fastq.gz LB2A_SRR1964642_trim.fastq.gz ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:45
Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
ILLUMINACLIP: Using 0 prefix pairs, 2 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Quality encoding detected as phred33
Input Reads: 26424138 Surviving: 25664909 (97.13%) Dropped: 759229 (2.87%)
TrimmomaticSE: Completed successfully



A total of 400mg of RNA was extracted from berry pericarp tissue (entire berries without seeds), a detailed description about RNA extraction and library preparation and sequencing of the 120 samples (10 varieties at four stages, in total 40 triplicate samples) can be also found in (Massonnet, M. et al. 2017)(https://pubmed.ncbi.nlm.nih.gov/28652263/). The complete information of the yielded 120 SRA ﬁles, downloaded from two BioProjects PRJNA265040 and PRJNA265039, can be found in supplemental tables. All the steps for RNA-seq analysis were performed in the UConn CBC Xanadu cluster belonging to the University of Connecticut, USA. This cluster also has SLURM as managing software. The general workﬂow used in this part can be observed in the workflow. For RNA data, after the quality control (QC) step, Trimmomatic was used to trim low quality and adapter contaminated sequences. In this case, the alignment of reads to the reference genome was performed by HISAT2 v2.2.1. HISAT2 is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. Before the alignment, the hisat2 build module was used to make a HISAT index ﬁle for the genome. By default, HISAT2 outputs the alignments in SAM format. Again samtools was used to sort the sequences, convert them to binary format and compress them. Finally, the function htseq-count from the HTSeq v0.13.5 package was used to count how many reads map to each annotated exon (gene) in the genome. The ﬁnal count for each gene was obtained from sum values for all their exons. These ﬁnal counts per gene are the inputs of the R package DESeq2 v3.13, used for the differential expression analysis.


## 2.3. FASTQC Before and After Quality Control
It is helpful to see how the quality of the data has changed after using Trimmomatic. To do this, we will be using the command-line versions of fastqc and MultiQC. These two programs create visual reports of the average quality of our reads.

dir="before"

module load fastqc/0.11.5
fastqc --outdir ./"$dir"/ ../raw_data/LB2A_SRR1964642.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LB2A_SRR1964643.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LC2A_SRR1964644.fastq.gz
fastqc --outdir ./"$dir"/ ../raw_data/LC2A_SRR1964645.fastq.gz
The full script for slurm scheduler is called fastqc_raw.sh and is located in the /fastqc folder.

The same command can be run on the fastq files after the trimming using fastqc program, and the comand will look like this:

dir="after"

module load fastqc/0.11.5
fastqc --outdir ./"$dir"/ ../quality_control/LB2A_SRR1964642.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/LB2A_SRR1964643.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/LC2A_SRR1964644.trim.fastq.gz -t 8
fastqc --outdir ./"$dir"/ ../quality_control/LC2A_SRR1964645.trim.fastq.gz -t 8
The full script for slurm scheduler is called fastqc_trimmed.sh which is located in /fastqc folder.

This will produce html files with the quality reports. The file strucutre inside the folder fastqc/ will look like this:

fastqc/
├── after
│   ├── trimmed_LB2A_SRR1964642_fastqc.html
│   ├── trimmed_LB2A_SRR1964642_fastqc.zip
│   ├── trimmed_LB2A_SRR1964643_fastqc.html
│   ├── trimmed_LB2A_SRR1964643_fastqc.zip
│   ├── trimmed_LC2A_SRR1964644_fastqc.html
│   ├── trimmed_LC2A_SRR1964644_fastqc.zip
│   ├── trimmed_LC2A_SRR1964645_fastqc.html
│   └── trimmed_LC2A_SRR1964645_fastqc.zip
├── before
│   ├── LB2A_SRR1964642_fastqc.html
│   ├── LB2A_SRR1964642_fastqc.zip
│   ├── LB2A_SRR1964643_fastqc.html
│   ├── LB2A_SRR1964643_fastqc.zip
│   ├── LC2A_SRR1964644_fastqc.html
│   ├── LC2A_SRR1964644_fastqc.zip
│   ├── LC2A_SRR1964645_fastqc.html
│   └── LC2A_SRR1964645_fastqc.zip
To view the html files you need to download them to your laptop and open them in a web browser. You can use a xanadu node dedicated to file transfer: transfer.cam.uchc.edu and the unix utility scp. Copy the files as shown below, or use an FTP client with a graphical user interface such as FileZilla or Cyberduck:

scp user-name@transfer.cam.uchc.edu:~/path/to/cloned/git/repository/fastqc/before/*.html .
The syntax is scp x y, meaning copy files x to location y. Do not forget the '.' at the end of the above code; which means to download the files to the current working directory in your computer. You can likewise download the HTML files for the trimmed reads.

Let's have a look at the output from fastqc. When loading the fastqc file, you will be greeted with this screen


There are some basic statistics which are all pretty self-explanatory. Notice that none of our sequence libraries fail the quality report! It would be concerning if we had even one because this report is from our trimmed sequence! The same thinking applies to our sequence length. Should the minimum of the sequence length be below 45, we would know that Trimmomatic had not run properly. Let's look at the next index in the file:



This screen is simply a box-and-whiskers plot of our quality scores per base pair. Note that there is a large variance and lower mean scores (but still about in our desired range) for base pairs 1-5 and that sequence quality declines toward the 3' end of the read.



This figure shows the distribution of mean read qualities. You can see we have a peak at about 38, which corresponds to a per base error probability of 0.00016.

The last panel at which we are going to look is the "Overrepresented Sequences" panel:


This is simply a list of sequences which appear disproportionately in our reads file. FastQC checks these against common adapter sequences and will flag them as such if they match. It is often the case in RNA-Seq that sequence from very highly expressed genes turns up in this panel. It may be helpful to try to identify some of these sequences using BLAST if they take up a large percentage of your library.

When you have a large experiment with many samples, checking FastQC HTML files can be a tedious task. To get around this, you can use use a program called MultiQC to combine them into a single report.

For HTML files in the before/ folder:

module load MultiQC/1.1




# 3. Variant Calling

# 4. eQTL analysis
