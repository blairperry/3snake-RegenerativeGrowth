# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# biocLite('IHW')
library(DESeq2)
library(IHW)

##
## Python
##
py.rawPSM <- read.csv('/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/PSM/Py_Raw/Python.MergedUnnormPSM_CleanID.csv',header=T,row.names = 'Accession')

x <- py.rawPSM
x <- x[,-1]

keep <- rowSums(x > 2) >= 2
x <- x[keep,]

hist(x[,2:8],nclass=20)

colData <- (DataFrame(condition=group <- factor(c('0hr','0hr','0hr','24hr','24hr','24hr','96hr','96hr'))))

dds <- DESeqDataSetFromMatrix(x,colData,formula(~ condition))
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)

Py_NormPSM <- counts(dds,normalized=TRUE)
#write.csv(as.data.frame(Py_NormPSM),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/PSM/Py.DEseq2NormPSM_11Apr18.csv')

#0v24
deRes0v24 <- as.data.frame(results(dds, contrast=c('condition','24hr','0hr')))
deRes0v24 <- deRes0v24[order(deRes0v24$padj),]

#write.csv(as.data.frame(deRes0v24),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Py_SI_0v24_DEseq2PSM_fil2_11Apr18.csv')


#24v96
deRes24v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','24hr')))
deRes24v96 <- deRes24v96[order(deRes24v96$padj),]

#write.csv(as.data.frame(deRes24v96),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Py_SI_24v96_DEseq2PSM_fil2_11Apr18.csv')


#0v96
deRes0v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','0hr')))
deRes0v96 <- deRes0v96[order(deRes0v96$padj),]

#write.csv(as.data.frame(deRes0v96),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Py_SI_0v96_DEseq2PSM_fil2_11Apr18.csv')

Py.Sig.0v24 <- sum(deRes0v24$padj < 0.1,na.rm=T)
Py.Sig.24v96 <- sum(deRes24v96$padj < 0.1,na.rm=T)
Py.Sig.0v96 <- sum(deRes0v96$padj < 0.1,na.rm=T)

##
## Watersnake
##

ner.rawPSM <- read.csv('/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/CV_Ner_RawTables_18Sep17/Ner/PSM/Ner.MergedUnnormPSM.csv',header=T,row.names = 'Accession')

x <- ner.rawPSM
x <- x[,-1]

keep <- rowSums(x > 2) >= 2
x <- x[keep,]

hist(x[,2:8],nclass=20)

colData <- (DataFrame(condition=group <- factor(c('0hr','0hr','0hr','24hr','24hr','24hr','96hr','96hr','96hr'))))

dds <- DESeqDataSetFromMatrix(x,colData,formula(~ condition))
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)

Ner_NormPSM <- counts(dds,normalized=TRUE)
#write.csv(as.data.frame(Ner_NormPSM),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/PSM/Ner.DEseq2NormPSM_11Apr18.csv')

#0v24
deRes0v24 <- as.data.frame(results(dds, contrast=c('condition','24hr','0hr')))
deRes0v24 <- deRes0v24[order(deRes0v24$padj),]

#write.csv(as.data.frame(deRes0v24),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Ner_SI_0v24_DEseq2PSM_fil2_11Apr18.csv')


#24v96
deRes24v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','24hr')))
deRes24v96 <- deRes24v96[order(deRes24v96$padj),]

#write.csv(as.data.frame(deRes24v96),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Ner_SI_24v96_DEseq2PSM_fil2_11Apr18.csv')


#0v96
deRes0v96 <- as.data.frame(results(dds, contrast=c('condition','96hr','0hr')))
deRes0v96 <- deRes0v96[order(deRes0v96$padj),]

#write.csv(as.data.frame(deRes0v96),file='/Users/perryb/Drive/_Research/_Projects/_3Species_Intestine_MS/data/proteomics/__NewProcessing/NormPSM_Stats/DEseq2_PairwiseResults/Ner_SI_0v96_DEseq2PSM_fil2_11Apr18.csv')

Ner.Sig.0v24 <- sum(deRes0v24$padj < 0.1,na.rm=T)
Ner.Sig.24v96 <- sum(deRes24v96$padj < 0.1,na.rm=T)
Ner.Sig.0v96 <- sum(deRes0v96$padj < 0.1,na.rm=T)



###
### Comparing numbers of DE proteins
###
library(viridis)
library(eulerr)
library(VennDiagram)
library(Vennerable)
library(latticeExtra)


DE.prots <- as.data.frame(c('fast.vs.1dpf','1dpf.v.4dpf','fast.v.4dpf'))
row.names(DE.prots) <- DE.prots[,1]

DE.prots$Python <- c(Py.Sig.0v24,Py.Sig.24v96,Py.Sig.0v96)
DE.prots$Watersnake <- c(Ner.Sig.0v24,Ner.Sig.24v96,Ner.Sig.0v96)
DE.prots <- DE.prots[,-1]

par(mfrow=c(1,2))
barplot(DE.prots[,1],main='Python',ylim=c(0,120),ylab='# of DE Proteins (Adj. p < 0.1', names.arg = c('Fasted vs. 1DPF', '1DPF vs. 4DPF','Fasted vs. 4DPF'),col='grey40')
barplot(DE.prots[,2],main='Watersnake',ylim=c(0,120), names.arg = c('Fasted vs. 1DPF', '1DPF vs. 4DPF','Fasted vs. 4DPF'),col='yellow1')
par(mfrow=c(1,1))



