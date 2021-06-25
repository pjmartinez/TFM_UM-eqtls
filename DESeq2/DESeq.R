
#####instalation of packages...just one time!

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("apeglm")
# BiocManager::install("pheatmap")
#BiocManager::install("vsn")
BiocManager::install("airway")
#############

# Load the libraries we'll need in the following code:
library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")

#White: Garganega       Glera   MoscatoBianco   Passerina       Vermentino
####17 de marzo 2021
directory <- "/Users/pedromartinez/Desktop/TFM_data/White_counts/soft_touch/" ###pathtofiles
sampleFiles_white <- list.files(directory, pattern = "*.counts_contados", full.names = F) 
View(sampleFiles_white)
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

sampleCondition_white <- c("Soft",
                     "Soft",
                     "Soft",
                     "Touch",
                     "Touch",
                     "Touch",
                     "Soft",
                     "Soft",
                     "Soft",
                     "Touch",
                     "Touch",
                     "Touch",
                     "Soft",
                     "Soft",
                     "Soft",
                     "Touch",
                     "Touch",
                     "Touch",
                     "Soft",
                     "Soft",
                     "Soft",
                     "Touch",
                     "Touch",
                     "Touch",
                     "Soft",
                     "Soft",
                     "Soft",
                     "Touch",
                     "Touch",
                     "Touch")



sampleTable_white <- data.frame(
  sampleName = sampleNames_white,
  fileName = sampleFiles_white,
  condition = sampleCondition_white
)



ddsHTSeq_white <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_white, 
  directory = directory, 
  design = ~ condition 
)


uva_white <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_white, 
  directory = directory, 
  design = ~ condition 
)


uva_white_test <- DESeq(uva_white)

resultadossinfiltrar_white  <- results(uva_white_test)

summary(resultadossinfiltrar_white)


sizeFactors(uva_white_test)

sumcounts_uva_white <- rowSums(counts(uva_white_test))


logsumcounts_uva_white <- log(sumcounts_uva_white,base=10)

hist(logsumcounts_uva_white,breaks=100)
keep_uva_white  <- sumcounts_uva_white  > 20
ddsHTSeq_uva_white_20filter <- uva_white_test[keep_uva_white,]

dds_uva_white_filter20 <- DESeq(ddsHTSeq_uva_white_20filter)
res_uva_white_filter20 <- results(dds_uva_white_filter20)
summary(res_uva_white_filter20)

normalized_counts_uva_white_filter20 <- counts(dds_uva_white_filter20, normalized=TRUE)

View(normalized_counts_uva_white_filter20)
write.table(normalized_counts_uva_white_filter20, file="normalized_FINAL_white_counts.txt", sep="\t", quote=F, col.names=NA)

summary(res_uva_white_filter20)
res_uva_white_filter20_05 <- results(dds_uva_white_filter20, alpha = 0.05)
summary(res_uva_white_filter20_05)

sum(res_uva_white_filter20_05$padj <0.05, na.rm = T)


res_shrink_res_uva_white_filter20  <- lfcShrink(dds_uva_white_filter20 ,coef="condition_Touch_vs_Soft")

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


uva_white_top20 <- order(-abs(res_shrink_res_uva_white_filter20$log2FoldChange))[1:20]
res_shrink_res_uva_white_filter20[uva_white_top20,]

plotMA(res_shrink_res_uva_white_filter20, ylim=c(-4,4))

# normalized, variance-stabilized transformed counts for visualization
vsd_dds_uva_white_filter20 <- vst(dds_uva_white_filter20, blind=FALSE)

plotPCA(vsd_dds_uva_white_filter20, intgroup="condition")





#######RED expression
####26 abril 2021
directory_red <- "/Users/pedromartinez/Desktop/TFM_data/Red_counts/touch_soft_red/counts/" ### estan en counts
sampleFiles_red <- list.files(directory_red, pattern = "*.counts_contados", full.names = F)
View(sampleFiles_red)
sampleNames_red <- c(
  "Barbera_Soft_1.counts_contados",
  "Barbera_Soft_2.counts_contados",
  "Barbera_Soft_3.counts_contados",
  "Barbera_Touch_1.counts_contados",
  "Barbera_Touch_2.counts_contados",
  "Barbera_Touch_3.counts_contados",
  "Negroamaro_Soft_1.counts_contados",
  "Negroamaro_Soft_2.counts_contados",
  "Negroamaro_Soft_3.counts_contados",
  "Negroamaro_Touch_1.counts_contados",
  "Negroamaro_Touch_2.counts_contados",
  "Negroamaro_Touch_3.counts_contados",
  "Primitivo_Soft_1.counts_contados",
  "Primitivo_Soft_2.counts_contados",
  "Primitivo_Soft_3.counts_contados",
  "Primitivo_Touch_1.counts_contados",
  "Primitivo_Touch_2.counts_contados",
  "Primitivo_Touch_3.counts_contados",
  "Refosco_Soft_1.counts_contados",
  "Refosco_Soft_2.counts_contados",
  "Refosco_Soft_3.counts_contados",
  "Refosco_Touch_1.counts_contados",
  "Refosco_Touch_2.counts_contados",
  "Refosco_Touch_3.counts_contados",
  "Sangiovese_Soft_1.counts_contados",
  "Sangiovese_Soft_2.counts_contados",
  "Sangiovese_Soft_3.counts_contados",
  "Sangiovese_Touch_1.counts_contados",
  "Sangiovese_Touch_2.counts_contados",
  "Sangiovese_Touch_3.counts_contados")

sampleCondition_red <- c("Soft",
                           "Soft",
                           "Soft",
                           "Touch",
                           "Touch",
                           "Touch",
                           "Soft",
                           "Soft",
                           "Soft",
                           "Touch",
                           "Touch",
                           "Touch",
                           "Soft",
                           "Soft",
                           "Soft",
                           "Touch",
                           "Touch",
                           "Touch",
                           "Soft",
                           "Soft",
                           "Soft",
                           "Touch",
                           "Touch",
                           "Touch",
                           "Soft",
                           "Soft",
                           "Soft",
                           "Touch",
                           "Touch",
                           "Touch")



sampleTable_red <- data.frame(
  sampleName = sampleNames_red,
  fileName = sampleFiles_red,
  condition = sampleCondition_red
)



uva_red <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_red, 
  directory = directory_red, 
  design = ~ condition 
)
uva_red$condition
treatments  <- c( "Soft","Touch")

uva_red$condition <- factor(uva_red$condition, levels=treatments)

uva_red_test <- DESeq(uva_red)
resultadossinfiltrar_red  <- results(uva_red_test)

summary(resultadossinfiltrar_red)

sizeFactors(uva_red_test)

sumcounts_uva_red <- rowSums(counts(uva_red_test))


logsumcounts_uva_red <- log(sumcounts_uva_red,base=10)

hist(logsumcounts_uva_red,breaks=100)
keep_uva_red  <- sumcounts_uva_red  > 20
uva_red_20filter <- uva_red_test[keep_uva_red,]

dds_uva_red_filter20 <- DESeq(uva_red_20filter)
res_uva_red_filter20 <- results(dds_uva_red_filter20)

s

summary(res_uva_red_filter20 )
res_uva_red_filter20_05 <- results(dds_uva_red_filter20, alpha = 0.05)
summary(res_uva_red_filter20_05)

sum(res_uva_red_filter20_05$padj <0.05, na.rm = T)

#results(uva_red_test)
Vitvi06g00229

normalized_counts_uva_red_filter20 <- counts(dds_uva_red_filter20, normalized=TRUE)

View(normalized_counts_uva_red_filter20)
write.table(normalized_counts_uva_red_filter20, file="normalized_FINAL_red_counts.txt", sep="\t", quote=F, col.names=NA)




res_shrink_dds_uva_red_filter20 <- lfcShrink(dds_uva_red_filter20 ,coef="condition_Touch_vs_Soft")

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=res_shrink_dds_uva_red_filter20$log2FoldChange,
  y=res_shrink_dds_uva_red_filter20$log2FoldChange,pch=20,
  cex=.2,
  col=1+(res_shrink_dds_uva_red_filter20$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change"
)
abline(0,1)


uva_red_top20 <- order(-abs(res_shrink_dds_uva_red_filter20$log2FoldChange))[1:20]
res_shrink_dds_uva_red_filter20[uva_red_top20,]

plotMA(res_shrink_dds_uva_red_filter20, ylim=c(-4,4))

# normalized, variance-stabilized transformed counts for visualization
vsd_dds_uva_red_filter20<- vst(dds_uva_red_filter20, blind=FALSE)

plotPCA(vsd_dds_uva_red_filter20, intgroup="condition")



#############28 abril ANALISIS GENERAL 
directory_analisis_general <- "/Users/pedromartinez/Desktop/TFM_data/analisis_general/" ### estan en counts
sampleFiles_analisis_general<- list.files(directory_analisis_general, pattern = "*.counts_contados", full.names = F)
View(sampleFiles_red)
sampleNames_analisis_general <- c(
  "Barbera_Soft_1.counts_contados",
  "Barbera_Soft_2.counts_contados",
  "Barbera_Soft_3.counts_contados",
  "Barbera_Touch_1.counts_contados",
  "Barbera_Touch_2.counts_contados",
  "Barbera_Touch_3.counts_contados",
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
  "Negroamaro_Soft_1.counts_contados",
  "Negroamaro_Soft_2.counts_contados",
  "Negroamaro_Soft_3.counts_contados",
  "Negroamaro_Touch_1.counts_contados",
  "Negroamaro_Touch_2.counts_contados",
  "Negroamaro_Touch_3.counts_contados",
  "Passerina_Soft_1.counts_contados",
  "Passerina_Soft_2.counts_contados",
  "Passerina_Soft_3.counts_contados",
  "Passerina_Touch_1.counts_contados",
  "Passerina_Touch_2.counts_contados",
  "Passerina_Touch_3.counts_contados",
  "Primitivo_Soft_1.counts_contados",
  "Primitivo_Soft_2.counts_contados",
  "Primitivo_Soft_3.counts_contados",
  "Primitivo_Touch_1.counts_contados",
  "Primitivo_Touch_2.counts_contados",
  "Primitivo_Touch_3.counts_contados",
  "Refosco_Soft_1.counts_contados",
  "Refosco_Soft_2.counts_contados",
  "Refosco_Soft_3.counts_contados",
  "Refosco_Touch_1.counts_contados",
  "Refosco_Touch_2.counts_contados",
  "Refosco_Touch_3.counts_contados",
  "Sangiovese_Soft_1.counts_contados",
  "Sangiovese_Soft_2.counts_contados",
  "Sangiovese_Soft_3.counts_contados",
  "Sangiovese_Touch_1.counts_contados",
  "Sangiovese_Touch_2.counts_contados",
  "Sangiovese_Touch_3.counts_contados",
  "Vermentino_Soft_1.counts_contados",
  "Vermentino_Soft_2.counts_contados",
  "Vermentino_Soft_3.counts_contados",
  "Vermentino_Touch_1.counts_contados",
  "Vermentino_Touch_2.counts_contados",
  "Vermentino_Touch_3.counts_contados")

sampleCondition_analisis_general <- c("Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch",
                         "Soft",
                         "Soft",
                         "Soft",
                         "Touch",
                         "Touch",
                         "Touch")



sampleTable_analisisgeneral <- data.frame(
    sampleName = sampleNames_analisis_general,
  fileName = sampleFiles_analisis_general,
  condition = sampleCondition_analisis_general
)



uva_analisis_general <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_analisisgeneral, 
  directory = directory_analisis_general, 
  design = ~ condition 
)


uva_analisis_general_test <- DESeq(uva_analisis_general)
resultadossinfiltrar_uvageneral  <- results(uva_analisis_general_test)

summary(resultadossinfiltrar_uvageneral)
?results

sizeFactors(uva_analisis_general_test)

sumcounts_analisis_general <- rowSums(counts(uva_analisis_general_test))


logsumcounts_uva_analisis_general  <- log(sumcounts_analisis_general,base=10)

hist(logsumcounts_uva_analisis_general,breaks=100)
keep_sumcounts_analisis_general  <- sumcounts_analisis_general  > 20
uva_analisis_general_20filter <- uva_analisis_general_test[keep_sumcounts_analisis_general,]

dds_uva_analisis_general_20filter <- DESeq(uva_analisis_general_20filter)
res_uva_analisis_general_20filter <- results(dds_uva_analisis_general_20filter)



summary(res_uva_analisis_general_20filter )
res_uva_analisis_general_20filter_05 <- results(dds_uva_analisis_general_20filter, alpha = 0.05)
summary(res_uva_analisis_general_20filter_05)

sum(res_uva_red_filter20_05$padj <0.05, na.rm = T)

#results(uva_red_test)
Vitvi06g00229

normalized_counts_uva_analisis_general_filter20 <- counts(dds_uva_analisis_general_20filter, normalized=TRUE)

View(normalized_counts_uva_analisis_general_filter20)
write.table(normalized_counts_uva_analisis_general_filter20, file="normalized_FINAL_ANALISISGENERAL_counts.txt", sep="\t", quote=F, col.names=NA)




res_shrink_dds_uva_analisis_general_20filter <- lfcShrink(dds_uva_analisis_general_20filter  ,coef="condition_Touch_vs_Soft")

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=res_shrink_dds_uva_analisis_general_20filter$log2FoldChange,
  y=res_shrink_dds_uva_analisis_general_20filter$log2FoldChange,pch=20,
  cex=.2,
  col=1+(res_shrink_dds_uva_analisis_general_20filter$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change"
)
abline(0,1)


uva_general_top20 <- order(-abs(res_shrink_dds_uva_analisis_general_20filter$log2FoldChange))[1:20]
res_shrink_dds_uva_analisis_general_20filter[uva_general_top20,]

plotMA(res_shrink_dds_uva_analisis_general_20filter, ylim=c(-4,4))




vsd_dds_uva_analisis_general_20filter<- vst(dds_uva_analisis_general_20filter, blind=FALSE)

plotPCA(vsd_dds_uva_analisis_general_20filter, intgroup="condition")



