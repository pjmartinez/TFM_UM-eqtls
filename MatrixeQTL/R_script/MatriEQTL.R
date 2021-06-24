install.packages("MatrixEQTL")
library(MatrixEQTL)
library(ggplot2)
setwd("/Users/pedromartinez/Desktop/TFM_data/")
basedir="/Users/pedromartinez/Desktop/TFM_data/TODOANTERIOR/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_wsoft = paste(basedir, "white_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_file_name_wsoft = paste(basedir, "mean_bysample_soft_white.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

errorCovariance = numeric();
errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps_wsoft = SlicedData$new();
snps_wsoft$fileDelimiter = "\t";      # the TAB character
snps_wsoft$fileOmitCharacters = "NA"; # denote missing values;
snps_wsoft$fileSkipRows = 1;          # one row of column labels
snps_wsoft$fileSkipColumns = 1;       # one column of row labels
snps_wsoft$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps_wsoft$LoadFile( SNP_file_name_wsoft );
summary(SNP_file_name_wsoft)

## Load gene expression data

gene_wsoft = SlicedData$new();
gene_wsoft$fileDelimiter = '\t'; # the TAB character
gene_wsoft$fileOmitCharacters = 'NA'; # denote missing values;
gene_wsoft$fileSkipRows = 1; # one row of column labels
gene_wsoft$fileSkipColumns = 1; # one column of row labels
gene_wsoft$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_wsoft$LoadFile(expression_file_name_wsoft);

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


#toc_load = proc.time()[3];
#Finally, the main Matrix eQTL function is called:



# ## Run the analysis
# {
#   tic_eqtl = proc.time()[3];
#   Matrix_eQTL_engine(snps, gene,cvrt, output_file_name,pvOutputThreshold,useModel, errorCovariance, verbose=TRUE,
#                      pvalue.hist = TRUE,
#                      min.pv.by.genesnp = FALSE,
#                      noFDRsaveMemory = FALSE);
#   toc_eqtl = proc.time()[3];
# }
# # cat('eQTL time: ', toc_eqtl-tic_eqtl, ' sec\n');
# cat('\n\n');
# show(data.frame(load = toc_load-tic_load, eQTL = toc_eqtl-tic_eqtl))
# 
#   
# me = Matrix_eQTL_engine(
#     snps = snps,
#     gene = gene,
#     cvrt = cvrt,
#     output_file_name = output_file_name,
#     pvOutputThreshold = pvOutputThreshold,
#     useModel = useModel,
#     errorCovariance = errorCovariance,
#     verbose = TRUE,
#     pvalue.hist = TRUE,
#     min.pv.by.genesnp = FALSE,
#     noFDRsaveMemory = FALSE);




######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis_wsoft = 'eQTL_cis_results_white_soft_ULTIMO_REPETICION_R.txt';
#output_file_name_tra = 'eQTL_tras_results_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
white_soft = Matrix_eQTL_main(
  snps = snps_wsoft,
  gene = gene_wsoft,
  cvrt = cvrt,
  useModel = useModel,
  pvOutputThreshold = 1e-9,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis_wsoft,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

red_soft$param



white_soft_REPETICION = Matrix_eQTL_main(
  snps = snps_wsoft,
  gene = gene_wsoft,
  cvrt = cvrt,
  useModel = useModel,
  pvOutputThreshold = 1e-9,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis_wsoft,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(me)
## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me$cis$eqtls$pvalue)
plot(white_soft, pch = 16, cex = 0.7, ylim = c(0,15));

dfFull = me$param$dfFull;
tstat = me$cis$eqtls$statistic;
r = tstat / sqrt( dfFull + tstat^2 );
R2 = r^2;

plot(me$cis$hist.counts)

me$all$eqtls$beta_se = me$all$eqtls$beta / me$all$eqtls$statistic;
NC_012017.3_860044_771558
View( white_soft$cis$eqtls)
Vitvi11g01337


cis_eqtl_res_white_soft = white_soft$cis$eqtls
cis_eqtl_res_white_soft = cis_eqtl_res_white_soft[cis_eqtl_res_white_soft$FDR < 0.1,]
top_eqtls_white_soft = cis_eqtl_res_white_soft[order(cis_eqtl_res_white_soft$pvalue),]
top_eqtls_white_soft = top_eqtls_white_soft[!duplicated(top_eqtls_white_soft$gene),]
top_eqtls_white_soft = merge(top_eqtls_white_soft, mafs, by="snps")
top_eqtls_white_soft = top_eqtls_white_soft[order(top_eqtls_white_soft$FDR),]
head(top_eqtls_white_soft,10)
View(top_eqtls_white_soft)
gene_id_white_soft = "Vitvi11g01337"
snp_id_white_soft =  "NC_012017.3_860044_771558"
#gene_id_white_soft = "Vitvi11g00497"
#snp_id_white_soft =  "NC_012017.3_5097457_896729"
snps_values_white= read.table("/Users/pedromartinez/Desktop/TFM_data/TODOANTERIOR/white_gt_1_0_2_noNA_nohomo.txt", row.names = 1, header = T)
gene_values= read.table("/Users/pedromartinez/Desktop/TFM_data/TODOANTERIOR/mean_bysample_soft_white.txt", row.names = 1, header = T)
snps_values_white=read.table("/Users/pedromartinez/Desktop/TFM_data/TODOANTERIOR/white_3_classes", row.names = 1, header = T)

# Get gene name of gene with lowest association p-value
#gene_id_white_soft = top_eqtls_white_soft$gene[6]
# Get corresponding SNP
snp_id_white_soft = top_eqtls_white_soft[top_eqtls_white_soft$gene == gene_id_white_soft,"snps"][1]
#snp_id_white_soft = cis_eqtl_res_white_soft[cis_eqtl_res_white_soft$gene == gene_id_white_soft,"snps"][2]
data_white_soft = data.frame(t(snps_values_white[snp_id_white_soft,]), t(gene_values[gene_id_white_soft,]))
# Get reference and alternative allele of the SNP
View(snps_values_white)
snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)#
View(snppos)
ref_alt_white_soft = unlist(snppos[snppos$CHR_POS_0 == snp_id_white_soft, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_white_soft= c(paste(ref_alt_white_soft[1], ref_alt_white_soft[1], sep="/"), paste(ref_alt_white_soft[1],
                                                                                         ref_alt_white_soft[2], sep="/"), paste(ref_alt_white_soft[2], ref_alt_white_soft[2], sep="/"))
gt_states_white_soft= factor(gt_states_white_soft, levels=gt_states_white_soft)
# Assign the labels
data_white_soft$gt = gt_states_white_soft[data_white_soft[,snp_id_white_soft]+1]
#data_white_touch$gt = gt_states_white_touch[data_white_touch[,snp_id_white_touch]+1]


# Subset to only genotype labels and expression
data_white_soft = data_white_soft[,c("gt", gene_id_white_soft)]
colnames(data_white_soft) = c("genotype", "expression")
# Plot
library(ggplot2)
white_soft_plot = ggplot(data_white_soft, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_white_soft, "with",snp_id_white_soft, "White Soft MODIFIER"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(white_soft_plot)






########################WHITE TOUCH CIS ANALYSIS

basedir="/Users/pedromartinez/Desktop/TFM_data/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_white = paste(basedir, "white_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_white_touch = paste(basedir, "mean_touch_white.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

errorCovariance = numeric();
errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps$LoadFile(SNP_file_name_white);


## Load gene expression data

gene_wtouch = SlicedData$new();
gene_wtouch$fileDelimiter = '\t'; # the TAB character
gene_wtouch$fileOmitCharacters = 'NA'; # denote missing values;
gene_wtouch$fileSkipRows = 1; # one row of column labels
gene_wtouch$fileSkipColumns = 1; # one column of row labels
gene_wtouch$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_wtouch$LoadFile(expression_white_touch);

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


#toc_load = proc.time()[3];
#Finally, the main Matrix eQTL function is called:



# ## Run the analysis
# {
#   tic_eqtl = proc.time()[3];
#   Matrix_eQTL_engine(snps, gene,cvrt, output_file_name,pvOutputThreshold,useModel, errorCovariance, verbose=TRUE,
#                      pvalue.hist = TRUE,
#                      min.pv.by.genesnp = FALSE,
#                      noFDRsaveMemory = FALSE);
#   toc_eqtl = proc.time()[3];
# }
# # cat('eQTL time: ', toc_eqtl-tic_eqtl, ' sec\n');
# cat('\n\n');
# show(data.frame(load = toc_load-tic_load, eQTL = toc_eqtl-tic_eqtl))
# 
#   
# me = Matrix_eQTL_engine(
#     snps = snps,
#     gene = gene,
#     cvrt = cvrt,
#     output_file_name = output_file_name,
#     pvOutputThreshold = pvOutputThreshold,
#     useModel = useModel,
#     errorCovariance = errorCovariance,
#     verbose = TRUE,
#     pvalue.hist = TRUE,
#     min.pv.by.genesnp = FALSE,
#     noFDRsaveMemory = FALSE);




######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis = 'eQTL_cis_results_white_touch_ULTIMOS_R.txt';
output_file_name_tra = 'eQTL_tras_results_white_touch_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
white_touch = Matrix_eQTL_main(
  snps = snps,
  gene = gene_wtouch,
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

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(white_touch)
## Results:

cat('Analysis done in: ', white_touch$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(white_touch$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)
View(white_touch$cis$eqtls)
## Plot the Q-Q plot of local and distant p-values

plot(white_touch$cis$eqtls$pvalue)
plot(white_touch, pch = 16, cex = 0.7, ylim = c(0,15))

###grabar sesion

save.image(file="eqtls_analisis_ultimos.RData")
load("/Users/pedromartinez/Desktop/TFM_data/eqtls_analisis_ultimos.RData")


######################


cis_eqtl_res_white_touch = white_touch$cis$eqtls
cis_eqtl_res_white_touch = cis_eqtl_res_white_touch[cis_eqtl_res_white_touch$FDR < 0.1,]
top_eqtls_white_touch = cis_eqtl_res_white_touch[order(cis_eqtl_res_white_touch$pvalue),]
top_eqtls_white_touch = top_eqtls_white_touch[!duplicated(top_eqtls_white_touch$gene),]
top_eqtls_white_touch = merge(top_eqtls_white_touch, mafs, by="snps")
top_eqtls_white_touch = top_eqtls_white_touch[order(top_eqtls_white_touch$FDR),]
head(top_eqtls_white_touch)



# Get gene name of gene with lowest association p-value
gene_id_white_touch = top_eqtls_white_touch$gene[3]
# Get corresponding SNP
snp_id_white_touch = top_eqtls_white_touch[top_eqtls_white_touch$gene == gene_id_white_touch,"snps"][1]
data_white_touch = data.frame(t(snps_values[snp_id_white_touch,]), t(gene_values[gene_id_white_touch,]))
# Get reference and alternative allele of the SNP

#snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)#
head(snppos)
ref_alt_white_touch = unlist(snppos[snppos$CHR_POS_0 == snp_id_white_touch, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_white_touch= c(paste(ref_alt_white_touch[1], ref_alt_white_touch[1], sep="/"), paste(ref_alt_white_touch[1],
                                                                                         ref_alt_white_touch[2], sep="/"), paste(ref_alt_white_touch[2], ref_alt_white_touch[2], sep="/"))
gt_states_white_touch = factor(gt_states_white_touch, levels=gt_states_white_touch)
# Assign the labels
data_white_touch$gt = gt_states_white_touch[data_white_touch[,snp_id_white_touch]+1]
# Subset to only genotype labels and expression
data_white_touch = data_white_touch[,c("gt", gene_id_white_touch)]
colnames(data_white_touch) = c("genotype", "expression")
# Plot
#library(ggplot2)
white_touch_plot = ggplot(data_white_touch, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_white_touch, "with",snp_id_white_touch, "White Touch"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(white_touch_plot)







########################RED TOUCH CIS ANALYSIS
set.seed(1234)
basedir="/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_red = paste(basedir, "red_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_red_touch = paste(basedir, "mean_bysample_red_touch.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

#errorCovariance = numeric();
#errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps_rtouch = SlicedData$new();
snps_rtouch$fileDelimiter = "\t";      # the TAB character
snps_rtouch$fileOmitCharacters = "NA"; # denote missing values;
snps_rtouch$fileSkipRows = 1;          # one row of column labels
snps_rtouch$fileSkipColumns = 1;       # one column of row labels
snps_rtouch$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps_rtouch$LoadFile(SNP_file_name_red);
#snps_wtouch

## Load gene expression data

gene_rtouch = SlicedData$new();
gene_rtouch$fileDelimiter = '\t'; # the TAB character
gene_rtouch$fileOmitCharacters = 'NA'; # denote missing values;
gene_rtouch$fileSkipRows = 1; # one row of column labels
gene_rtouch$fileSkipColumns = 1; # one column of row labels
gene_rtouch$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_rtouch$LoadFile(expression_red_touch);

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



######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis = 'eQTL_cis_results_red_touch_ULTIMO_R.txt';
#output_file_name_tra = 'eQTL_tras_results_red_touch_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
red_touch = Matrix_eQTL_main(
  snps = snps_rtouch,
  gene = gene_rtouch,
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

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(red_touch)
## Results:

cat('Analysis done in: ', white_touch$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(red_touch$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(red_touch$cis$eqtls$pvalue)
plot(red_touch, pch = 16, cex = 0.7, ylim = c(0,15))

View(red_touch$cis$eqtls)


cis_eqtl_res_red_touch = red_touch$cis$eqtls
cis_eqtl_res_red_touch = cis_eqtl_res_red_touch[cis_eqtl_res_red_touch$FDR < 0.1,]
top_eqtls_red_touch = cis_eqtl_res_red_touch[order(cis_eqtl_res_red_touch$pvalue),]
top_eqtls_red_touch = top_eqtls_red_touch[!duplicated(top_eqtls_red_touch$gene),]
top_eqtls_red_touch = merge(top_eqtls_red_touch, mafs, by="snps")
top_eqtls_red_touch = top_eqtls_red_touch[order(top_eqtls_red_touch$FDR),]
head(top_eqtls_red_touch)

NC_012017.3_860044_771558
gene_id_white_soft = "Vitvi11g01337"
NC_012021.3_12757455_5405454

# Get gene name of gene with lowest association p-value
gene_id_red_touch = top_eqtls_red_touch$gene[5]
Vitvi15g01453
gene_id_red_touch ="Vitvi15g01453"
gene_id_red_touch ="Vitvi05g00925"
# Get corresponding SNP
snp_id_red_touch = top_eqtls_red_touch[top_eqtls_red_touch$gene == gene_id_red_touch,"snps"][1]

snp_id_red_touch="NC_012021.3_12757455_5405454"
data_red_touch = data.frame(t(snps_values_redsoft[snp_id_red_touch,]), t(gene_values[gene_id_red_touch,]))
# Get reference and alternative allele of the SNP

#snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)#
head(snppos)
ref_alt_red_touch = unlist(snppos[snppos$CHR_POS_0 == snp_id_red_touch, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_red_touch= c(paste(ref_alt_red_touch[1], ref_alt_red_touch[1], sep="/"), paste(ref_alt_red_touch[1],
                                                                                      ref_alt_red_touch[2], sep="/"), paste(ref_alt_red_touch[2], ref_alt_red_touch[2], sep="/"))
gt_states_red_touch = factor(gt_states_red_touch, levels=gt_states_red_touch)
# Assign the labels
data_red_touch$gt = gt_states_red_touch[data_red_touch[,snp_id_red_touch]+1]
# Subset to only genotype labels and expression
data_red_touch = data_red_touch[,c("gt", gene_id_red_touch)]
colnames(data_red_touch) = c("genotype", "expression")
# Plot
#library(ggplot2)
red_touch_plot = ggplot(data_red_touch, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_red_touch, "with",snp_id_red_touch, "Red Touch"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(red_touch_plot)








save.image(file="eqtls_analisis_ultimos.RData")
load("/Users/pedromartinez/Desktop/TFM_data/eqtls_analisis_ultimos.RData")




#######################RED SOFT CIS ANALYSIS

basedir="/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_red = paste(basedir, "red_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_red_soft = paste(basedir, "mean_bysample_red_soft.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

#errorCovariance = numeric();
#errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps_rsoft = SlicedData$new();
snps_rsoft$fileDelimiter = "\t";      # the TAB character
snps_rsoft$fileOmitCharacters = "NA"; # denote missing values;
snps_rsoft$fileSkipRows = 1;          # one row of column labels
snps_rsoft$fileSkipColumns = 1;       # one column of row labels
snps_rsoft$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps_rsoft$LoadFile(SNP_file_name_red);


## Load gene expression data

gene_rsoft = SlicedData$new();
gene_rsoft$fileDelimiter = '\t'; # the TAB character
gene_rsoft$fileOmitCharacters = 'NA'; # denote missing values;
gene_rsoft$fileSkipRows = 1; # one row of column labels
gene_rsoft$fileSkipColumns = 1; # one column of row labels
gene_rsoft$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_rsoft$LoadFile(expression_red_soft);

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




######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis = 'eQTL_cis_results_red_soft_ULTIMO_R.txt';
#output_file_name_tra = 'eQTL_tras_results_red_soft_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
red_soft = Matrix_eQTL_main(
  snps = snps_rsoft,
  gene = gene_rsoft,
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

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(red_soft)
## Results:

cat('Analysis done in: ', red_soft$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(red_soft$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(red_soft$cis$eqtls$pvalue)
plot(red_soft, pch = 16, cex = 0.7, ylim = c(0,15))


geneshigh <- top_eqtls_red_soft[top_eqtls_red_soft$gene =="Vitvi06g00964",]


"/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/mean_bysample_red_soft.tx"
gene_values_redsoft= read.table("/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/mean_bysample_red_soft.txt", row.names = 1, header = T)
snps_values_redsoft= read.table("/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/red_gt_1_0_2_noNA_nohomo.txt", row.names = 1, header = T)

View(geneshigh)

cis_eqtl_res_red_soft = red_soft$cis$eqtls
cis_eqtl_res_red_soft = cis_eqtl_res_red_soft[cis_eqtl_res_red_soft$FDR < 0.1,]
top_eqtls_red_soft = cis_eqtl_res_red_soft[order(cis_eqtl_res_red_soft$pvalue),]
top_eqtls_red_soft = top_eqtls_red_soft[!duplicated(top_eqtls_red_soft$gene),]
top_eqtls_red_soft = merge(top_eqtls_red_soft, mafs, by="snps")
top_eqtls_red_soft = top_eqtls_red_soft[order(top_eqtls_red_soft$FDR),]
head(top_eqtls_red_soft)
12942414
Vitvi06g00964

top_eqtls_red_soft$snps ==
VIT_08s0040g00860
Vitvi08g02114

# Get gene name of gene with lowest association p-value
gene_id_red_soft = top_eqtls_red_soft$gene[4]
gene_id_red_soft = "Vitvi11g00497"
gene_id_red_soft = "Vitvi06g00964"
snp_id_red_soft = "NC_012012.3_12942414_14288906"
# Get corresponding SNP
snp_id_red_soft = top_eqtls_red_soft[top_eqtls_red_soft$gene == gene_id_red_soft,"snps"][1]
data_red_soft = data.frame(t(snps_values_redsoft[snp_id_red_soft,]), t(gene_values_redsoft[gene_id_red_soft,]))
# Get reference and alternative allele of the SNP

#snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)#
head(snppos)
ref_alt_red_soft = unlist(snppos[snppos$CHR_POS_0 == snp_id_red_soft, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_red_soft= c(paste(ref_alt_red_soft[1], ref_alt_red_soft[1], sep="/"), paste(ref_alt_red_soft[1],
                                                                                                  ref_alt_red_soft[2], sep="/"), paste(ref_alt_red_soft[2], ref_alt_red_soft[2], sep="/"))
gt_states_red_soft = factor(gt_states_red_soft, levels=gt_states_red_soft)
# Assign the labels
data_red_soft$gt = gt_states_red_soft[data_red_soft[,snp_id_red_soft]+1]
# Subset to only genotype labels and expression
data_red_soft = data_red_soft[,c("gt", gene_id_red_soft)]
colnames(data_red_soft) = c("genotype", "expression")
# Plot
#library(ggplot2)
red_soft_plot = ggplot(data_red_soft, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_red_soft, "with",snp_id_red_soft, "Red Soft", "HIGH IMPACT"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(red_soft_plot)










save.image(file="eqtls_analisis_ultimos.RData")

#######################GENERAL ANALYSIS SOFT

basedir="/Users/pedromartinez/Desktop/TFM_data/analisis_general/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_General = paste(basedir, "reorder_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_soft_General = paste(basedir, "mean_total_soft.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

#errorCovariance = numeric();
#errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

snps_sfgeneral = SlicedData$new();
snps_sfgeneral$fileDelimiter = "\t";      # the TAB character
snps_sfgeneral$fileOmitCharacters = "NA"; # denote missing values;
snps_sfgeneral$fileSkipRows = 1;          # one row of column labels
snps_sfgeneral$fileSkipColumns = 1;       # one column of row labels
snps_sfgeneral$fileSliceSize = 10000;      # read file in pieces of 2,000 rows
snps_sfgeneral$LoadFile(SNP_file_name_General);


## Load gene expression data

gene_sfgeneral = SlicedData$new();
gene_sfgeneral$fileDelimiter = '\t'; # the TAB character
gene_sfgeneral$fileOmitCharacters = 'NA'; # denote missing values;
gene_sfgeneral$fileSkipRows = 1; # one row of column labels
gene_sfgeneral$fileSkipColumns = 1; # one column of row labels
gene_sfgeneral$fileSliceSize = 10000; # read file in pieces of 10,000 rows
gene_sfgeneral$LoadFile(expression_soft_General);

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




######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis = 'eQTL_cis_results_GENERAL_soft_ULTIMO_R.txt';
#output_file_name_tra = 'eQTL_tras_results_GENERAL_soft_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
GENERAL_soft = Matrix_eQTL_main(
  snps = snps_sfgeneral,
  gene = gene_sfgeneral,
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

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(red_soft)
## Results:

cat('Analysis done in: ', red_soft$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(GENERAL_soft$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(GENERAL_soft$cis$eqtls$pvalue)
plot(GENERAL_soft, pch = 16, cex = 0.7, ylim = c(0,25))





cis_eqtl_res_GENERAL_soft = GENERAL_soft$cis$eqtls
cis_eqtl_res_GENERAL_soft = cis_eqtl_res_GENERAL_soft[cis_eqtl_res_GENERAL_soft$FDR < 0.1,]
top_eqtls_GENERAL_soft = cis_eqtl_res_GENERAL_soft[order(cis_eqtl_res_GENERAL_soft$pvalue),]
top_eqtls_GENERAL_soft = top_eqtls_GENERAL_soft[!duplicated(top_eqtls_GENERAL_soft$gene),]
mafs = apply(as.matrix(snps_values),1,mean)/2
mafs = data.frame(snps=names(mafs), maf = mafs)
top_eqtls_GENERAL_soft = merge(top_eqtls_GENERAL_soft, mafs, by="snps")
top_eqtls_GENERAL_soft = top_eqtls_GENERAL_soft[order(top_eqtls_GENERAL_soft$FDR),]
head(top_eqtls_GENERAL_soft)

NC_012020.3_3829243_3761383

# Get gene name of gene with lowest association p-value
gene_id_GENERAL_soft = top_eqtls_GENERAL_soft$gene[1]
snp_id_GENERAL_soft="NC_012023.3_4085526_6909284"
gene_id_GENERAL_soft = "Vitvi17g00347"


gene_id_GENERAL_soft = "Vitvi14g00297"


gene_id_GENERAL_soft = "Vitvi09g02012"
snp_id_GENERAL_soft="NC_012015.3_22447315_17258222"

gene_id_GENERAL_soft = "Vitvi17g00347"
snp_id_GENERAL_soft="NC_012023.3_4085526_6909284"

gene_id_GENERAL_soft = "Vitvi01g01498"
snp_id_GENERAL_soft="NC_012007.3_20337582_10371804"



gene_id_GENERAL_soft = "Vitvi14g00297"
snp_id_GENERAL_soft="NC_012020.3_3829243_3761383"
gene_id_GENERAL_soft = "Vitvi17g00347"
snp_id_GENERAL_soft="NC_012023.3_4085526_6909284"
gene_id_GENERAL_soft = "Vitvi01g01498"
snp_id_GENERAL_soft="NC_012007.3_20337582_10371804"
gene_id_GENERAL_soft = "Vitvi03g01373"
snp_id_GENERAL_soft="NC_012009.3_1043545_11131276"
gene_id_GENERAL_soft = "Vitvi05g01272"
snp_id_GENERAL_soft="NC_012011.3_18911113_13576026"
gene_id_GENERAL_soft = "Vitvi09g01367"
snp_id_GENERAL_soft="NC_012015.3_21267175_17202317"
gene_id_GENERAL_soft = "Vitvi09g02012"
snp_id_GENERAL_soft="NC_012015.3_22440993_17257909"

Vitvi14g00297
# Get corresponding SNP
#snp_id_GENERAL_soft = top_eqtls_GENERAL_soft[top_eqtls_GENERAL_soft$gene == gene_id_GENERAL_soft,"snps"][1]
data_GENERAL_soft = data.frame(t(snps_values[snp_id_GENERAL_soft,]), t(gene_values[gene_id_GENERAL_soft,]))
# Get reference and alternative allele of the SNP
#snp_id_GENERAL_soft="NC_012020.3_3829243_3761383"
#snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)
head(snppos)
ref_alt_GENERAL_soft = unlist(snppos[snppos$CHR_POS_0 == snp_id_GENERAL_soft, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_GENERAL_soft= c(paste(ref_alt_GENERAL_soft[1], ref_alt_GENERAL_soft[1], sep="/"), paste(ref_alt_GENERAL_soft[1],
                                                                                                     ref_alt_GENERAL_soft[2], sep="/"), paste(ref_alt_GENERAL_soft[2], ref_alt_GENERAL_soft[2], sep="/"))
gt_states_GENERAL_soft = factor(gt_states_GENERAL_soft, levels=gt_states_GENERAL_soft)
# Assign the labels
data_GENERAL_soft$gt = gt_states_GENERAL_soft[data_GENERAL_soft[,snp_id_GENERAL_soft]+1]
# Subset to only genotype labels and expression
data_GENERAL_soft = data_GENERAL_soft[,c("gt", gene_id_GENERAL_soft)]
colnames(data_GENERAL_soft) = c("genotype", "expression")
# Plot

soft = ggplot(data_GENERAL_soft, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_GENERAL_soft, "with",snp_id_GENERAL_soft, "Soft General p-value=3.82e-13"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(soft)







#######################GENERAL ANALYSIS TOUCH

basedir="/Users/pedromartinez/Desktop/TFM_data/analisis_general/"

useModel = modelLINEAR; # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name_General = paste(basedir, "reorder_gt_1_0_2_noNA_nohomo.txt", sep="");
expression_touch_General = paste(basedir, "mean_total_touch.txt", sep="");

covariates_file_name = paste(basedir, "Covariates2.txt", sep="");
covariates_file_name = character();

#The p-value threshold determines which gene-SNP associations are saved in the output file output_file_name. Note that for larger datasets the threshold should be lower. Setting the threshold to a high value for a large dataset may cause excessively large output files.

#pvOutputThreshold = 1e-2;

#Finally, define the covariance matrix for the error term. This parameter is rarely used. If the covariance matrix is a multiple of identity, set it to numeric().

#errorCovariance = numeric();
#errorCovariance = character();
#The next section of the sample code contains three very similar parts loading the files with genotype, gene expression, and covariates. In each part one can set the file delimiter (i.e. tabulation "\t", comma ",", or space " "), the string representation for missing values, the number of rows with column labels, and the number of columns with row labels. Finally, one can change the number of the variables in a slice for the file reading procedure (do not change if not sure).

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
gene_tchgeneral$LoadFile(expression_touch_General);

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




######CIS TRAS


snps_location_file_name = paste(basedir, "snpsloc.txt", sep="");
gene_location_file_name = paste(basedir, "geneloc.txt", sep="");

# Output file name
#Output file name
output_file_name_cis = 'eQTL_cis_results_GENERAL_touch_ULTIMO_R.txt';
#output_file_name_tra = 'eQTL_tras_results_GENERAL_touch_R.txt';

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-8;
#pvOutputThreshold_tra = 1e-8;

# Distance for local gene-SNP pairs
cisDist = 1e3;


## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

set.seed(1234)
GENERAL_touch = Matrix_eQTL_main(
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

#unlink(output_file_name_tra);
#unlink(output_file_name_cis);

str(GENERAL_touch)
## Results:

cat('Analysis done in: ', GENERAL_touch$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(GENERAL_touch$cis$eqtls)
#cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(GENERAL_touch$cis$eqtls$pvalue)
plot(GENERAL_touch, pch = 16, cex = 0.7, ylim = c(0,25))


#snps_values= read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", row.names = 1, header = T)
gene_values= read.table("/Users/pedromartinez/Desktop/TFM_data/analisis_general/mean_total_touch.txt", row.names = 1, header = T)
snps_values= read.table("/Users/pedromartinez/Desktop/TFM_data/analisis_general/reorder_gt_1_0_2_noNA_nohomo.txt", row.names = 1, header = T)

cis_eqtl_res_GENERAL_touch = GENERAL_touch$cis$eqtls
cis_eqtl_res_GENERAL_touch = cis_eqtl_res_GENERAL_touch[cis_eqtl_res_GENERAL_touch$FDR < 0.1,]
top_eqtls_GENERAL_touch = cis_eqtl_res_GENERAL_touch[order(cis_eqtl_res_GENERAL_touch$pvalue),]
top_eqtls_GENERAL_touch = top_eqtls_GENERAL_touch[!duplicated(top_eqtls_GENERAL_touch$gene),]
mafs = apply(as.matrix(snps_values),1,mean)/2
mafs = data.frame(snps=names(mafs), maf = mafs)
top_eqtls_GENERAL_touch = merge(top_eqtls_GENERAL_touch, mafs, by="snps")
top_eqtls_GENERAL_touch = top_eqtls_GENERAL_touch[order(top_eqtls_GENERAL_touch$FDR),]
head(top_eqtls_GENERAL_touch)

VIT_07s0005g06300
Vitvi07g00905

gene_id_GENERAL_touch = "Vitvi12g02656"
snp_id_GENERAL_touch = "NC_012018.3_19505811_2465034"
gene_id_GENERAL_touch = "Vitvi15g01618"
snp_id_GENERAL_touch = "NC_012021.3_18081894_5616606"
gene_id_GENERAL_touch = "Vitvi19g01396"
snp_id_GENERAL_touch = "NC_012025.3_17470744_9324677"
gene_id_GENERAL_touch = "Vitvi01g01520"
snp_id_GENERAL_touch = "NC_012007.3_20595239_10381474"
gene_id_GENERAL_touch = "Vitvi04g01283"
snp_id_GENERAL_touch = "NC_012010.3_18571715_12667623"
gene_id_GENERAL_touch = "Vitvi05g01272"
snp_id_GENERAL_touch = "NC_012011.3_18911113_13576026"
gene_id_GENERAL_touch = "Vitvi08g00849"
snp_id_GENERAL_touch = "NC_012014.3_10631405_15825489"

# Get gene name of gene with lowest association p-value
gene_id_GENERAL_touch = top_eqtls_GENERAL_touch$gene[84]
gene_id_GENERAL_touch = "Vitvi12g02656"
Vitvi12g02656

#snp_id_GENERAL_touch = "NC_012018.3_19505811_2465034"
#gene_id_GENERAL_touch = "Vitvi01g01552"

#gene_id_GENERAL_touch = "Vitvi02g00189"
# Get corresponding SNP
#snp_id_GENERAL_touch = top_eqtls_GENERAL_touch[top_eqtls_GENERAL_touch$gene == gene_id_GENERAL_touch,"snps"][1]
data_GENERAL_touch = data.frame(t(snps_values[snp_id_GENERAL_touch,]), t(gene_values[gene_id_GENERAL_touch,]))
# Get reference and alternative allele of the SNP

#snppos <- read.table("/Users/pedromartinez/Desktop/TFM_data/plink_TFM_data/res_snps_refall", sep="\t", header=TRUE)
head(snppos)
ref_alt_GENERAL_touch = unlist(snppos[snppos$CHR_POS_0 == snp_id_GENERAL_touch, c("REF", "ALT")])
# Prepare the genotype labels
gt_states_GENERAL_touch= c(paste(ref_alt_GENERAL_touch[1], ref_alt_GENERAL_touch[1], sep="/"), paste(ref_alt_GENERAL_touch[1],
                                                           ref_alt_GENERAL_touch[2], sep="/"), paste(ref_alt_GENERAL_touch[2], ref_alt_GENERAL_touch[2], sep="/"))
gt_states_GENERAL_touch = factor(gt_states_GENERAL_touch, levels=gt_states_GENERAL_touch)
# Assign the labels
data_GENERAL_touch$gt = gt_states_GENERAL_touch[data_GENERAL_touch[,snp_id_GENERAL_touch]+1]
# Subset to only genotype labels and expression
data_GENERAL_touch = data_GENERAL_touch[,c("gt", gene_id_GENERAL_touch)]
colnames(data_GENERAL_touch) = c("genotype", "expression")
# Plot
library(ggplot2)
p = ggplot(data_GENERAL_touch, aes(genotype, expression)) +
  ggtitle(paste("eQTL of gene",gene_id_GENERAL_touch, "with",snp_id_GENERAL_touch, "Touch General p-value=6.69e-11"))+
  geom_jitter(colour="darkgrey", position=position_jitter(width=0.25)) +
  geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + theme_bw()
print(p)

