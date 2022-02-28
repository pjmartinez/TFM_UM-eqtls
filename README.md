# TFM_UM-eqtls
## Title: Algorithms for the discovery of cis-eQTL signals in woody species: the vine (*Vitis vinifera* L.) as a study model.


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a Slurm scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
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

## 2. Accessing the Data using SRA-Toolkit
Before we can get started, we need to get the data we're going to analyze. This dataset has been deposited in the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra) at NCBI, a comprehensive collection of sequenced genetic data submitted by researchers. The beauty of the SRA is the ease with which genetic data becomes accessible to any scientist with an internet connection. Sets of sequences (usually all the sequences from a given sample within an experiment) in the SRA have a unique identifier. The set may be downloaded using a software module called the `sratoolkit`. There are a variety of commands in the `sratoolkit`, which I invite you to investigate for yourself at [here](https://www.ncbi.nlm.nih.gov/books/NBK569238/).

An overview of the project data can be viewed here.

We will download data from the control samples (LB2A) and the heat stress treatment (LC2A). The SRA accessions are as follows:

LB2A : SRR1964642, SRR1964643
LC2A : SRR1964644, SRR1964645

We have provided a script to download the data from the the SRA data using SRA-Toolkit using this script. It contains three sections.

The SLURM header:

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

module load sratoolkit/2.8.1 

fastq-dump --gzip SRR1964642
fastq-dump --gzip SRR1964643
fastq-dump --gzip SRR1964644
fastq-dump --gzip SRR1964645
The line module load sratoolkit/2.8.1 loads the sratoolkit so we can use it. Once downloaded, we rename the samples:

mv SRR1964642.fastq.gz LB2A_SRR1964642.fastq.gz
mv SRR1964643.fastq.gz LB2A_SRR1964643.fastq.gz
mv SRR1964644.fastq.gz LC2A_SRR1964644.fastq.gz
mv SRR1964645.fastq.gz LC2A_SRR1964645.fastq.gz
The full script for slurm scheduler can be found in the raw_data folder. Before running it, add your own e-mail address to the --mail-user option (or delete the line entirely if you don't want an e-mail notification when the job completes).

When you're ready, you can execute the script by entering sbatch fastq_dump_xanadu.sh in the terminal. This submits the job to the SLURM scheduler.

Once the job is completed the folder structure will look like this:

raw_data/
|-- fastq_dump_xanadu_NNNN.err
|-- fastq_dump_xanadu_NNNN.out
|-- LB2A_SRR1964642.fastq.gz
|-- LB2A_SRR1964643.fastq.gz
|-- LC2A_SRR1964644.fastq.gz
`-- LC2A_SRR1964645.fastq.gz
The sequence files are in fastq format and compressed using gzip (indicated by the .gz). It's good practice to keep sequence files compressed. Most bioinformatics programs can read them directly without needing to decompress them first, and it doesn't get in the way of inspecting them either. The .out and .err files are output produced by SLURM that you can use to troubleshoot if things go wrong. Lets have a look at at the contents of one of the fastq files:

zcat LB2A_SRR1964642.fastq.gz | head -n 12

@SRR1964642.1 FCC355RACXX:2:1101:1476:2162 length=90
CAACATCTCAGTAGAAGGCGGCGCCTTCACCTTCGACGTGGGGAATCGCTTCAACCTCACGGGGGCTTTCCTCTACACGTCCTGTCCGGA
+SRR1964642.1 FCC355RACXX:2:1101:1476:2162 length=90
?@@D?DDBFHHFFGIFBBAFG:DGHDFHGHIIIIC=D<:?BBCCCCCBB@BBCCCB?CCBB<@BCCCAACCCCC>>@?@88?BCACCBB>
@SRR1964642.2 FCC355RACXX:2:1101:1641:2127 length=90
NGCCTGTAAAATCAAGGCATCCCCTCTCTTCATGCACCTCCTGAAATAAAAGGGCCTGAATAATGTCGTACAGAAGACTGCGGCACAGAC
+SRR1964642.2 FCC355RACXX:2:1101:1641:2127 length=90
#1=DDFFFHHHHGJJJJJIIIJIJGIIJJJIJIJJGIJIJJJJIJJJJJJIJJJIJJJJJJJGIIHIGGHHHHHFFFFFDEDBDBDDDDD
@SRR1964642.3 FCC355RACXX:2:1101:1505:2188 length=90
GGACAACGCCTGGACTCTGGTTGGTATTGTCTCCTGGGGAAGCAGCCGTTGCTCCACCTCCACTCCTGGTGTCTATGCCCGTGTCACCGA
+SRR1964642.3 FCC355RACXX:2:1101:1505:2188 length=90
CCCFFFFFHHFFHJJJIIIJHHJJHHJJIJIIIJEHJIJDIJJIIJJIGIIIIJGHHHHFFFFFEEEEECDDDDEDEDDDDDDDADDDDD
Each sequence record has four lines. The first is the sequence name, beginning with @. The second is the nucleotide sequence. The third is a comment line, beginning with +, and which here only contains the sequence name again (it is often empty). The fourth are phred-scaled base quality scores, encoded by ASCII characters. Follow the links to learn more, but in short, the quality scores give the probability a called base is incorrect.

3. Quality Control Using Trimmomatic



A total of 400mg of RNA was extracted from berry pericarp tissue (entire berries without seeds), a detailed description about RNA extraction and library preparation and sequencing of the 120 samples (10 varieties at four stages, in total 40 triplicate samples) can be also found in (Massonnet, M. et al. 2017)(https://pubmed.ncbi.nlm.nih.gov/28652263/). The complete information of the yielded 120 SRA ﬁles, downloaded from two BioProjects PRJNA265040 and PRJNA265039, can be found in supplemental tables. All the steps for RNA-seq analysis were performed in the UConn CBC Xanadu cluster belonging to the University of Connecticut, USA. This cluster also has SLURM as managing software. The general workﬂow used in this part can be observed in the workflow. For RNA data, after the quality control (QC) step, Trimmomatic was used to trim low quality and adapter contaminated sequences. In this case, the alignment of reads to the reference genome was performed by HISAT2 v2.2.1. HISAT2 is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. Before the alignment, the hisat2 build module was used to make a HISAT index ﬁle for the genome. By default, HISAT2 outputs the alignments in SAM format. Again samtools was used to sort the sequences, convert them to binary format and compress them. Finally, the function htseq-count from the HTSeq v0.13.5 package was used to count how many reads map to each annotated exon (gene) in the genome. The ﬁnal count for each gene was obtained from sum values for all their exons. These ﬁnal counts per gene are the inputs of the R package DESeq2 v3.13, used for the differential expression analysis.


# 3. Variant Calling

# 4. eQTL analysis
