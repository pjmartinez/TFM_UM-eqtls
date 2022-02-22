# TFM_UM-eqtls
## Title: Algorithms for the discovery of cis-eQTL signals in woody species: the vine (*Vitis vinifera* L.) as a study model.


This repository is a publicly available tutorial for eQTL analysis using RNA-Seq data and DNA data in woody species. All steps should be run in a cluster with appropriate headers for a Slurm scheduler that can be modified simply to run. Commands should never be executed on the submit nodes of any HPC machine. Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs. If you are new to Linux, please use this handy guide for the operating system commands. In this guide, you will be working with common bio Informatic file formats, such as FASTA, FASTQ, SAM/BAM, and GFF3/GTF. You can learn even more about each file format here. 
In summary, the repository includs the different inputs, scripts, outputs files (the supplemental files) obtained during the analysis performance in the biorxv paper: https://www.biorxiv.org/content/10.1101/2021.07.06.450811v1





![Screenshot](/Figures/generalpipeline.png)



# 1. Overview
In this study, expresion transcriptional responses to three environmental stresses were examined in the large yellow croaker (Larimichthys crocea): heat stress (LC2A), cold stress (LA2A) and a 21-day fast (LF1A). mRNA was extracted from livers from from fishes from the three treatments and a control (LB2A) and sequenced on an Illumina HiSeq 2000. Sequencing was single-end, and each sequence is 90bp.

For the purposes of this tutorial, we will compare the control (LB2A) to the heat stress treatment (LC2A).

## Cloning the tutorial repository

To work through this tutorial, copy it to your home directory using the git clone command:
```
git clone < git-repository-path.git > 

```

