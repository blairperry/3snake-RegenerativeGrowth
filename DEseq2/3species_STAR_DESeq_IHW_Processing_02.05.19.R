#BiocManager::install("apeglm", version = "3.8")

library(dplyr)
library(DESeq2)
library(IHW)
library(ggplot2)
library(viridis)
library(ggrepel)
library(apeglm)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/STAR_Feb2019")

nonZeroGenes <- read.csv('3species_HumanIDNonZeroList_02.05.19.csv',header = T)

###
### PYTHON
###
py.rawCounts <- read.table('raw_counts/python_rawCounts_HumanIDs_02.05.19.txt',row.names=1,stringsAsFactors = F,sep='\t',header=F)
py.rawCounts <- py.rawCounts[,c(-1,-2,-3,-4)]

count(rowSums(py.rawCounts[,3:20]) == 0)

py.GeneToProt <- as.data.frame(row.names(py.rawCounts))
py.GeneToProt$protID <- py.rawCounts$protein_id

py.SnakeToHuman <- as.data.frame(row.names(py.rawCounts))
py.SnakeToHuman$humID <- py.rawCounts$V26

py.rawCounts.simple <- py.rawCounts[,c(-1,-2,-21)]

colNames <- as.data.frame(colnames(py.rawCounts.simple))

colData <- (DataFrame(condition=group <- factor(c('0', '0', '0', '0', '0', '0', '1', '1', '1', '1', '4', '4', '4', '4', '4', '4', '4', '4')),
                      type=group <- factor(c('a','b','a','a','b','b','a','b','b','a','a','b','b','a','b','a','b','a'))))


dds <- DESeqDataSetFromMatrix(py.rawCounts.simple,colData,formula(~ type + condition))

dds <- dds[rowSums(counts(dds)) > 1,]

dds <- DESeq(dds)
 
py_normcounts <- counts(dds,normalized=TRUE)

#write.csv(as.data.frame(py_normcounts),file='norm_counts/py_SI_normcounts_DEseq_02.05.2019.csv')

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="condition_1_vs_0", type="apeglm")


 
py_deRes0v24 <- results(dds, contrast=c('condition','1','0'))
plotMA(py_deRes0v24,ylim=c(-2,2))

plotMA(resLFC,ylim=c(-2,2))


py_deRes24v96 <- results(dds, contrast=c('condition','4','1'))
plotMA(py_deRes24v96)


vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition",'type'), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type, label=name)) +
  geom_point(size=3) +
  #geom_text() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()



##########
# 0v24
##########

deRes0v24 <- as.data.frame(results(dds, contrast=c('condition','1','0')))

ihwRes <- ihw(pvalue ~ baseMean,  data = deRes0v24, alpha = 0.05)

rejections(ihwRes)
deRes <- deRes0v24
deRes <- na.omit(deRes)



plot(ihwRes)

gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg
gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)

ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 10)

ggplot(deRes, aes(x=pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)

ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")

rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
                 covariate_type="log2 fc")) %>%
  ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
  facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

deRes0v24$IHW_pvalue <- ihwRes@df$adj_pvalue

count(deRes0v24$IHW_pvalue < 0.05,na.rm = T)

deRes0v24.wProtID <- merge(deRes0v24,py.GeneToProt,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wProtID <- deRes0v24.wProtID[order(deRes0v24.wProtID$IHW_pvalue),]

deRes0v24.wHumID <- merge(deRes0v24,py.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wHumID <- deRes0v24.wHumID[order(deRes0v24.wHumID$IHW_pvalue),]
deRes0v24.wHumID <- deRes0v24.wHumID[which(deRes0v24.wHumID$humID %in% nonZeroGenes[,1]),]

count(deRes0v24.wHumID$IHW_pvalue < 0.05,na.rm = T)

write.csv(as.data.frame(deRes0v24.wHumID),file='./pairwise_results/python_0v24_DEseq2Out_NonZeroFiltered_02.05.19.csv',row.names = F)


#024v96
deRes24v96 <- as.data.frame(results(dds, contrast=c('condition','4','1')))
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes24v96, alpha = 0.05)
rejections(ihwRes)
deRes <- deRes24v96
deRes <- na.omit(deRes)

plot(ihwRes)

gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
  geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
gg
gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)

ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)

deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 10)

ggplot(deRes, aes(x=pvalue)) +
  geom_histogram(binwidth = 0.025, boundary = 0) +
  facet_wrap( ~ baseMeanGroup, nrow = 2)
ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")

rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
                 covariate_type="base mean"),
      data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
                 covariate_type="log2 fc")) %>%
  ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
  facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

deRes24v96$IHW_pvalue <- ihwRes@df$adj_pvalue
deRes24v96.wProtID <- merge(deRes24v96,py.GeneToProt,by.x='row.names',by.y=1,all.x=T)

deRes24v96.wProtID <- deRes24v96.wProtID[order(deRes24v96.wProtID$IHW_pvalue),]

deRes24v96.wHumID <- merge(deRes24v96,py.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes24v96.wHumID <- deRes24v96.wHumID[order(deRes24v96.wHumID$IHW_pvalue),]
deRes24v96.wHumID <- deRes24v96.wHumID[which(deRes24v96.wHumID$humID %in% nonZeroGenes[,1]),]

write.csv(as.data.frame(deRes24v96.wHumID),file='./pairwise_results/python_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv',row.names = F)



##########################################################################################################################################################

###
### RATTLESNAKE
###
cv.rawCounts <- read.table('raw_counts/CVV_allSamples_rawCounts_HumanIDs_02.05.19.txt',row.names=1,stringsAsFactors = F,sep='\t',header=F)
cv.rawCounts <- cv.rawCounts[,c(-1,-2,-3,-4)]

count(rowSums(cv.rawCounts[,3:13]) == 0)


cv.GeneToProt <- as.data.frame(row.names(cv.rawCounts))
cv.GeneToProt$TransID <- cv.rawCounts$Crovir_Transcript_ID

cv.SnakeToHuman <- as.data.frame(row.names(cv.rawCounts))
cv.SnakeToHuman$humID <- cv.rawCounts$V19

cv.rawCounts.simple <- cv.rawCounts[,c(-1,-2,-14)]

colData <- (DataFrame(condition=group <- factor(c('0hr','0hr','0hr','0hr','24hr','24hr','24hr','96hr','96hr','96hr','96hr'))))


dds <- DESeqDataSetFromMatrix(cv.rawCounts.simple,colData,formula(~ condition))

dds <- dds[rowSums(counts(dds)) > 1,]

dds <- DESeq(dds)


vsd <- vst(dds, blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=3) +
  #geom_text() +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()




cv_normcounts <- counts(dds,normalized=TRUE)
#write.csv(as.data.frame(cv_normcounts),file='norm_counts/cv_SI_normcounts_DEseq_02.05.2019.csv')

R_deRes0v24 <- results(dds, contrast=c('condition','24hr','0hr'))
plotMA(R_deRes0v24)

R_deRes24v96 <- results(dds, contrast=c('condition','96hr','24hr'))
plotMA(R_deRes24v96)


#0v24
deRes0v24 <- as.data.frame(results(dds, contrast=c('condition','24hr','0hr')))
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes0v24, alpha = 0.05)

rejections(ihwRes)
deRes <- deRes0v24
deRes <- na.omit(deRes)

 # plot(ihwRes)
# 
# gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
#   geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
# gg
# gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)
# 
# ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
# 
# deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 14)
# 
# ggplot(deRes, aes(x=pvalue)) +
#   geom_histogram(binwidth = 0.025, boundary = 0) +
#   facet_wrap( ~ baseMeanGroup, nrow = 2)
# ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")
# 
# rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
#                  covariate_type="base mean"),
#       data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
#                  covariate_type="log2 fc")) %>%
#   ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
#   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))


deRes0v24$IHW_pvalue <- ihwRes@df$adj_pvalue
deRes0v24.wProtID <- merge(deRes0v24,cv.GeneToProt,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wProtID <- deRes0v24.wProtID[order(deRes0v24.wProtID$IHW_pvalue),]

deRes0v24.wHumID <- merge(deRes0v24,cv.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wHumID <- deRes0v24.wHumID[order(deRes0v24.wHumID$IHW_pvalue),]
deRes0v24.wHumID <- deRes0v24.wHumID[which(deRes0v24.wHumID$humID %in% nonZeroGenes[,1]),]

write.csv(as.data.frame(deRes0v24.wHumID),file='./pairwise_results/rattlesnake_0v24_DEseq2Out_NonzeroFiltered_02.05.19.csv',row.names = F)




#24v96
deRes24v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','24hr')))

ihwRes <- ihw(pvalue ~ baseMean,  data = deRes24v96, alpha = 0.05)

rejections(ihwRes)
deRes <- deRes24v96
deRes <- na.omit(deRes)

# plot(ihwRes)
# 
# gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
#   geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
# gg
# gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)
# 
# ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
# 
# deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 14)
# 
# ggplot(deRes, aes(x=pvalue)) +
#   geom_histogram(binwidth = 0.025, boundary = 0) +
#   facet_wrap( ~ baseMeanGroup, nrow = 2)
# 
# ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")
# 
# rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
#                  covariate_type="base mean"),
#       data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
#                  covariate_type="log2 fc")) %>%
#   ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
#   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))


deRes24v96$IHW_pvalue <- ihwRes@df$adj_pvalue
deRes24v96.wProtID <- merge(deRes24v96,cv.GeneToProt,by.x='row.names',by.y=1,all.x=T)
deRes24v96.wProtID <- deRes24v96.wProtID[order(deRes24v96.wProtID$IHW_pvalue),]

deRes24v96.wHumID <- merge(deRes24v96,cv.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes24v96.wHumID <- deRes24v96.wHumID[order(deRes24v96.wHumID$IHW_pvalue),]
deRes24v96.wHumID <- deRes24v96.wHumID[which(deRes24v96.wHumID$humID %in% nonZeroGenes[,1]),]

write.csv(as.data.frame(deRes24v96.wHumID),file='./pairwise_results/rattlesnake_24v96_DEseq2Out_NonzeroFiltered_02.05.19.csv',row.names = F)



##########################################################################################################################################################

###
### Watersnake
###
ner.rawCounts <- read.table('raw_counts/Ner_allSamples_rawCounts_HumanIDs_02.05.19.txt',row.names=1,stringsAsFactors = F,sep='\t')
ner.rawCounts <- ner.rawCounts[,c(-1,-2,-3,-4)]

count(rowSums(ner.rawCounts[,3:14]) == 0)

ner.GeneToProt <- as.data.frame(row.names(ner.rawCounts))
ner.GeneToProt$protID <- ner.rawCounts$V7

ner.SnakeToHuman <- as.data.frame(row.names(ner.rawCounts))
ner.SnakeToHuman$humID <- ner.rawCounts$V20

ner.rawCounts.simple <- ner.rawCounts[,c(-1,-2,-15)] 

colData <- (DataFrame(condition=group <- factor(c('0hr','0hr','0hr','0hr','24hr','24hr','24hr','24hr','96hr','96hr','96hr','96hr'))))

dds <- DESeqDataSetFromMatrix(ner.rawCounts.simple,colData,formula(~ condition))

dds <- dds[rowSums(counts(dds)) > 1,]

dds <- DESeq(dds)
 
ner_normcounts <- counts(dds,normalized=TRUE)

#write.csv(as.data.frame(ner_normcounts),file='norm_counts/ner_SI_normcounts_DEseq_02.05.2019.csv')

W_deRes0v24 <- results(dds, contrast=c('condition','24hr','0hr'))
plotMA(W_deRes0v24)
W_deRes24v96 <- results(dds, contrast=c('condition','96hr','24hr'))
plotMA(W_deRes24v96)

 
#0v24
deRes0v24 <- as.data.frame(results(dds, contrast=c('condition','24hr','0hr')))
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes0v24, alpha = 0.05)

rejections(ihwRes)
deRes <- deRes0v24
deRes <- na.omit(deRes)

# plot(ihwRes)
# 
# gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
#   geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
# gg
# gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)
# 
# ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
# 
# deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 13)
# 
# ggplot(deRes, aes(x=pvalue)) +
#   geom_histogram(binwidth = 0.025, boundary = 0) +
#   facet_wrap( ~ baseMeanGroup, nrow = 2)
# 
# ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")
# 
# rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
#                  covariate_type="base mean"),
#       data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
#                  covariate_type="log2 fc")) %>%
#   ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
#   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

deRes0v24$IHW_pvalue <- ihwRes@df$adj_pvalue
deRes0v24.wProtID <- merge(deRes0v24,ner.GeneToProt,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wProtID <- deRes0v24.wProtID[order(deRes0v24.wProtID$IHW_pvalue),]

deRes0v24.wHumID <- merge(deRes0v24,ner.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes0v24.wHumID <- deRes0v24.wHumID[order(deRes0v24.wHumID$IHW_pvalue),]
deRes0v24.wHumID <- deRes0v24.wHumID[which(deRes0v24.wHumID$humID %in% nonZeroGenes[,1]),]

write.csv(as.data.frame(deRes0v24.wHumID),file='./pairwise_results/watersnake_0v24_DEseq2Out_NonzeroFiltered_02.05.19.csv',row.names = F)




#024v96
deRes24v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','24hr')))
ihwRes <- ihw(pvalue ~ baseMean,  data = deRes24v96, alpha = 0.05)
rejections(ihwRes)

deRes <- deRes24v96
deRes <- na.omit(deRes)


# plot(ihwRes)
# 
# gg <- ggplot(as.data.frame(ihwRes), aes(x = pvalue, y = adj_pvalue, col = group)) +
#   geom_point(size = 0.25) + scale_colour_hue(l = 70, c = 150, drop = FALSE)
# gg
# gg %+% subset(as.data.frame(ihwRes), adj_pvalue <= 0.2)
# 
# ggplot(deRes, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
# 
# deRes$baseMeanGroup <- groups_by_filter(deRes$baseMean, 13)
# 
# ggplot(deRes, aes(x=pvalue)) +
#   geom_histogram(binwidth = 0.025, boundary = 0) +
#   facet_wrap( ~ baseMeanGroup, nrow = 2)
# 
# ggplot(deRes, aes(x = pvalue, col = baseMeanGroup)) + stat_ecdf(geom = "step")
# 
# rbind(data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$baseMean)/nrow(deRes),
#                  covariate_type="base mean"),
#       data.frame(pvalue = deRes$pvalue, covariate = rank(deRes$log2FoldChange)/nrow(deRes),
#                  covariate_type="log2 fc")) %>%
#   ggplot(aes(x = covariate, y = -log10(pvalue))) + geom_hex(bins = 100) +
#   facet_grid( . ~ covariate_type) + ylab(expression(-log[10]~p))

deRes24v96$IHW_pvalue <- ihwRes@df$adj_pvalue
deRes24v96.wProtID <- merge(deRes24v96,ner.GeneToProt,by.x='row.names',by.y=1,all.x=T)
deRes24v96.wProtID <- deRes24v96.wProtID[order(deRes24v96.wProtID$IHW_pvalue),]

deRes24v96.wHumID <- merge(deRes24v96,ner.SnakeToHuman,by.x='row.names',by.y=1,all.x=T)
deRes24v96.wHumID <- deRes24v96.wHumID[order(deRes24v96.wHumID$IHW_pvalue),]
deRes24v96.wHumID <- deRes24v96.wHumID[which(deRes24v96.wHumID$humID %in% nonZeroGenes[,1]),]

write.csv(as.data.frame(deRes24v96.wHumID),file='./pairwise_results/watersnake_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv',row.names = F)
