# TFM_UM-eqtls
## Title: Algorithms for the discovery of cis-eQTL signals in woody species: the vine (Vitis vinifera) as a study model.


This respository represents the complete pipeline to perfomance eQTLs analysis in grapes using MatrixeQTL. In addition, this repository host the different output files obtained during the analysis.


![Screenshot](/Figures/generalpipeline.png)


## Prepare a reference genome
The first step is the prepare the reference genome. Most software packages that align short-read sequencing data to, or otherwise manipulate a reference genome require that genome to be indexed in some way. We will generate indexes using both `bwa` and `samtools`. For the workshop, a pre-indexed human genome is provided, but the code below shows how you can download and index a human genome yourself:
