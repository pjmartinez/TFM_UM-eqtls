
#####instalacion de los paquetes para su uso...una sola vez

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
# create an object with the directory containing your counts:
# !!edit this to point to your own count file directory!!
# directory <- "/Users/pedromartinez/Desktop/TFM_data/comparacion_DEG/"
# 
# # ensure the count files are where you think they are
# list.files(directory)
# 
# sampleFiles <- list.files(directory, pattern = "*.counts", full.names = F)
# View(sampleFiles)
# 
# # sampleNames <- c(
#   "Barbera_Soft_1.counts",
#   "Barbera_Soft_2.counts",
#   "Barbera_Soft_3.counts",
#   "Barbera_Touch_1.counts",
#   "Barbera_Touch_2.counts",
#   "Barbera_Touch_3.counts",
#   "Garganega_Soft_1.counts",
#   "Garganega_Soft_2.counts",
#   "Garganega_Soft_3.counts",
#   "Garganega_Touch_1.counts",
#   "Garganega_Touch_2.counts",
#   "Garganega_Touch_3.counts",
#   "Glera_Soft_1.counts",
#   "Glera_Soft_2.counts",
#   "Glera_Soft_3.counts",
#   "Glera_Touch_1.counts",
#   "Glera_Touch_2.counts",
#   "Glera_Touch_3.counts",
#   "Moscatobianco_Soft_1.counts",
#   "Moscatobianco_Soft_2.counts",
#   "Moscatobianco_Soft_3.counts",
#   "Moscatobianco_Touch_1.counts",
#   "Moscatobianco_Touch_2.counts",
#   "Moscatobianco_Touch_3.counts",
#   "Negroamaro_Soft_1.counts",
#   "Negroamaro_Soft_2.counts",
#   "Negroamaro_Soft_3.counts",
#   "Negroamaro_Touch_1.counts",
#   "Negroamaro_Touch_2.counts",
#   "Negroamaro_Touch_3.counts",
#   "Passerina_Soft_1.counts",
#   "Passerina_Soft_2.counts",
#   "Passerina_Soft_3.counts",
#   "Passerina_Touch_1.counts",
#   "Passerina_Touch_2.counts",
#   "Passerina_Touch_3.counts",
#   "Primitivo_Soft_1.counts",
#   "Primitivo_Soft_2.counts",
#   "Primitivo_Soft_3.counts",
#   "Primitivo_Touch_1.counts",
#   "Primitivo_Touch_2.counts",
#   "Primitivo_Touch_3.counts",
#   "Refosco_Soft_1.counts",
#   "Refosco_Soft_2.counts",
#   "Refosco_Soft_3.counts",
#   "Refosco_Touch_1.counts",
#   "Refosco_Touch_2.counts",
#   "Refosco_Touch_3.counts",
#   "Sangiovese_Soft_1.counts",
#   "Sangiovese_Soft_2.counts",
#   "Sangiovese_Soft_3.counts",
#   "Sangiovese_Touch_1.counts",
#   "Sangiovese_Touch_2.counts",
#   "Sangiovese_Touch_3.counts",
#   "Vermentino_Soft_1.counts",
#   "Vermentino_Soft_2.counts",
#   "Vermentino_Soft_3.counts",
#   "Vermentino_Touch_1.counts",
#   "Vermentino_Touch_2.counts",
#   "Vermentino_Touch_3.counts")

#White: Garganega       Glera   MoscatoBianco   Passerina       Vermentino
# SNP_id 	1
# "../align_stepwise/Garganega_mkdup.bam"	4
# "../align_stepwise/Glera_mkdup.bam"	5
# "../align_stepwise/MoscatoBianco_mkdup.bam"	6
# "../align_stepwise/Passerina_mkdup.bam"	8
# "../align_stepwise/Vermentino_mkdup.bam"	11
####17 de marzo 2021
directory <- "/Users/pedromartinez/Desktop/TFM_data/White_counts/soft_touch/"
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





#######red expression
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



####CALCULATE MEAN VALUES


# ddsHTSeq2 <- DESeqDataSetFromHTSeqCount(
#   sampleTable = sampleTable, 
#   directory = directory, 
#   design = ~ condition + type
# )

# dds2 <- DESeq(ddsHTSeq2)

ddsHTSeq_white <- estimateSizeFactors(ddsHTSeq_white)

sizeFactors(ddsHTSeq_white)






normalized_counts_white2_test <- counts(ddsHTSeq_white, normalized=TRUE)
write.table(normalized_counts_white2_test, file="normalized_TEST2_white_counts.txt", sep="\t", quote=F, col.names=NA)

dds_white_test <- DESeq(ddsHTSeq_white)


View(dds_white)
ddsHTSeq_white$condition
results_white <- results(dds_white)
counts(dds_white)

# colData(ddsHTSeq2)$type<-factor(colData(ddsHTSeq2)$type, levels=c("White", "Red"))
# dds2 <- DESeq(ddsHTSeq2)
# res2 <- results(dds2)
# res
normalized_counts_white <- counts(dds_white, normalized=TRUE)

View(normalized_counts_white)
write.table(normalized_counts_white, file="normalized_TEST_white_counts.txt", sep="\t", quote=F, col.names=NA)

colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("Soft", "Touch"))

sumcounts_white <- rowSums(counts(ddsHTSeq_white))


logsumcounts_white <- log(sumcounts_white,base=10)

hist(logsumcounts_white,breaks=100)
keep_white <- sumcounts_white > 20
ddsHTSeq_white_20filter <- ddsHTSeq_white[keep_white,]

dds_white_filter20 <- DESeq(ddsHTSeq_white_20filter)
res_white_filter20 <- results(dds_white_filter20)


normalized_counts_white_filter20 <- counts(dds_white_filter20, normalized=TRUE)

View(normalized_counts_white_filter20)
write.table(normalized_counts_white_filter20, file="normalized_TEST_white_counts.txt", sep="\t", quote=F, col.names=NA)




View(res_white_filter20)



res<-res[order(res$padj),]
mcols(res,use.names=TRUE)
summary(res)
resultsNames(dds)
summary(res)
head(res)
dim(res)


#GENERATE MA PLOT
png('MAPlot_htseq.png')
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()

#WRITE RESULTS INTO FILE
write.csv(as.data.frame(res),file="deseq2_htseq_C1_vs_C2.csv")


plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

d<- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
              returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))



rld_white_filter20 <- rlog(dds_white_filter20, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld), 3)

ntd_white_filter20 <- normTransform(dds_white_filter20)
library("vsn")
notAllZero <- (rowSums(counts(dds_white_filter20))>0)
meanSdPlot(assay(ntd_white_filter20)[notAllZero,])


library("pheatmap")
select_white_filter20 <- order(rowMeans(counts(dds_white_filter20,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df_white_filter20 <- as.data.frame(colData(dds_white_filter20))
pheatmap(assay(ntd_white_filter20)[select_white_filter20,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_white_filter20)

pheatmap(assay(rld_white_filter20)[select_white_filter20,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df_white_filter20)

BiocManager::install("pasilla")
colData(dds)

library("pasilla")
data("pasillaExons")
data("pasillaGenes")
data("dxd")
browseVignettes("pasilla")


par(mar=c(14,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, cex.axis=0.5)

boxplot(assays(dds)[["cooks"]])


summary(res)
# get shrunken log fold changes
res_shrink_white_filter20 <- lfcShrink(dds_white_filter20,coef="condition_Touch_vs_Soft")

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=res_white_filter20$log2FoldChange,
  y=res_shrink_white_filter20$log2FoldChange,pch=20,
  cex=.2,
  col=1+(res_white_filter20$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change"
)
abline(0,1)

# get the top 20 genes by shrunken log2 fold change
top20 <- order(-abs(res_shrink$log2FoldChange))[1:20]
res_shrink[top20,]

plotMA(res_shrink, ylim=c(-4,4))



# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.2,
     col=(log_padj > 10)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

# normalized, variance-stabilized transformed counts for visualization
vsd_white_filter20 <- vst(dds_white_filter20, blind=FALSE)

plotPCA(vsd_white_filter20, intgroup="condition")

# alternatively, using ggplot

dat_white_filter20 <- plotPCA(vsd_white_filter20, intgroup="condition",returnData=TRUE)



View(dat$name)

type <- c("Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "Red",
          "White",
          "White",
          "White",
          "White",
          "White",
          "White")
dat$type <- type
dim(dat)
p <- ggplot(dat,aes(x=PC1,y=PC2,col=group, shape=type))
p + geom_point()




# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# get top 50 log fold change genes
top50 <- order(-abs(res_shrink$log2FoldChange))[1:50]
df <- data.frame(colData(dds)[,"condition"])
rownames(df) <- colnames(dds)
colnames(df) <- "condition"
pheatmap(
  assay(rld)[top50,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
)


l2fc_ord <- order(-abs(res_shrink$log2FoldChange))
plotCounts(dds, gene=l2fc_ord[1], intgroup="condition")




rld2 <- rlog(dds, blind=TRUE)
plotPCA(rld2, intgroup="condition")
rld_mat <- assay(rld2)
pca <- prcomp(t(rld_mat))


df2 <- as.data.frame(pca$x)
ggplot(df2) + geom_point(aes(x=PC3, y=PC4))


library("airway")
data("airway")
View(airway)
se <- airway
View(se)
