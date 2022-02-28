# TFM_UM-eqtls
## Title: Algorithms for the discovery of cis-eQTL signals in woody species: the vine (*Vitis vinifera* L.) as a study model.


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a Slurm scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
In summary, the repository includs the different inputs, scripts, outputs files (the supplemental files) obtained during the analysis performance in the biorxv paper: https://www.biorxiv.org/content/10.1101/2021.07.06.450811v1





![Screenshot](/Figures/generalpipeline.png)


Contest
 1. [Overview](README.md/Overwiew)
 2. [RNA Seq]
 3. [Variant Calling]
 4. [eQTL analysis]

# 1. Overview
This tutorial will teach you how to use open source quality control, RNA Seq, Variant Calling, eQTL tools to complete a cis-eQTL analysis which is possible when you  have generated the specific datasets. Moving through the tutorial, you will take expression and genotypic data from a woody species as grape and perform a eQTL analysis via Matrix eQTL to characterize the gene expression levels during ripening Vitis vinifera L. fruit.



1. PartA

2. PartB

3. PartC
    
the aim of this study is to integrate the gene expression information with genome structural features characterize the gene expression levels during ripening Vitis vinifera L. fruit. For that public DNA and RNA data from 10 grape cultivars, 5 red and 5 white cultivars, were used. The integration of both layers of  information, structural variants and differential expression, around ripening onset allowed a cis-eQTL analysis to identify genes and mechanisms relevant for this complex process. through the use of bioinformatics approaches commonly used in human research. The ﬁnal inventory of cis-eQTLs will be an important resource for future research to understand the mechanism for variation in gene regulation during ripening in this species, and could be considered general markers of ripening in grapes.


#RNA Seq Analysis

A total of 400mg of RNA was extracted from berry pericarp tissue (entire berries without seeds), a detailed description about RNA extraction and library preparation and sequencing of the 120 samples (10 varieties at four stages, in total 40 triplicate samples) can be also found in (Massonnet, M. et al. 2017). The complete information of the yielded 120 SRA ﬁles, downloaded from two BioProjects PRJNA265040 and PRJNA265039, can be found in supplemental tables. All the steps for RNA-seq analysis were performed in the UConn CBC Xanadu cluster belonging to the University of Connecticut, USA. This cluster also has SLURM as managing software. The general workﬂow used in this part can be observed in the workflow. For RNA data, after the quality control (QC) step, Trimmomatic was used to trim low quality and adapter contaminated sequences. In this case, the alignment of reads to the reference genome was performed by HISAT2 v2.2.1. HISAT2 is a fast and sensitive aligner for mapping next generation sequencing reads against a reference genome. Before the alignment, the hisat2 build module was used to make a HISAT index ﬁle for the genome. By default, HISAT2 outputs the alignments in SAM format. Again samtools was used to sort the sequences, convert them to binary format and compress them. Finally, the function htseq-count from the HTSeq v0.13.5 package was used to count how many reads map to each annotated exon (gene) in the genome. The ﬁnal count for each gene was obtained from sum values for all their exons. These ﬁnal counts per gene are the inputs of the R package DESeq2 v3.13, used for the differential expression analysis.

## Cloning the tutorial repository

To work through this tutorial, copy it to your home directory using the git clone command:
```
git clone < git-repository-path.git > 

```

