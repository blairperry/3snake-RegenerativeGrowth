# install.packages('eulerr')
# install.packages('viridis')
# install.packages('VennDiagram')
# source("https://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph",'reshape','gtools','xtable'))
# library(devtools)
# BiocManager::install("RBGL", version = "3.8")
# BiocManager::install("graph", version = "3.8")
# install_github("js229/Vennerable")
# install.packages('latticeExtra')

library(eulerr)
library(viridis)
library(VennDiagram)
library(Vennerable)
library(latticeExtra)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/")

### Pairwise Significant Genes - 0v24 p-val < 0.05

py_0v24data = read.csv("STAR_Feb2019/pairwise_results/python_0v24_DEseq2Out_NonZeroFiltered_02.05.19.csv", row.names=1)
cv_0v24data = read.csv("STAR_Feb2019/pairwise_results/rattlesnake_0v24_DEseq2Out_NonzeroFiltered_02.05.19.csv", row.names=1)
ner_0v24data = read.csv("STAR_Feb2019/pairwise_results/watersnake_0v24_DEseq2Out_NonzeroFiltered_02.05.19.csv", row.names=1)

py_0v24data = py_0v24data[which(py_0v24data$IHW_pvalue < 0.05),]
cv_0v24data = cv_0v24data[which(cv_0v24data$IHW_pvalue < 0.05),]
ner_0v24data = ner_0v24data[which(ner_0v24data$IHW_pvalue < 0.05),]

py_0v24sigIDs = py_0v24data$humID 
cv_0v24sigIDs = cv_0v24data$humID 
ner_0v24sigIDs = ner_0v24data$humID

x = list('py'=py_0v24sigIDs,'cv'=cv_0v24sigIDs,'ner'=ner_0v24sigIDs)

vennD=Venn(x)
vennD

#Enter second row of vennD result, starting at value below '100'
fit0v24 <- euler(c('A'=1601,'B'=670,'A&B'=563,'C'=296,'A&C'=191,'B&C'=102,'A&B&C'=204))

Venn0v24 = plot(fit0v24,main='0v24 Pairwise Significant Genes',
                fill=viridis(3),
                fill_opacity = .5,
                auto.key=T,
                quantities=T,
                labels=c('Python','Rattlesnake','Watersnake'))

Venn0v24

### Pairwise Significant Genes - 24v96 p-val < 0.05

py_24v96data = read.csv("STAR_Feb2019/pairwise_results/python_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv", row.names=1)
cv_24v96data = read.csv("STAR_Feb2019/pairwise_results/rattlesnake_24v96_DEseq2Out_NonzeroFiltered_02.05.19.csv", row.names=1)
ner_24v96data = read.csv("STAR_Feb2019/pairwise_results/watersnake_24v96_DEseq2Out_NonzeroFiltered_02.05.19.csv", row.names=1)

py_24v96data = py_24v96data[which(py_24v96data$IHW_pvalue < 0.05),]
cv_24v96data = cv_24v96data[which(cv_24v96data$IHW_pvalue < 0.05),]
ner_24v96data = ner_24v96data[which(ner_24v96data$IHW_pvalue < 0.05),]

py_24v96sigIDs = py_24v96data$humID 
cv_24v96sigIDs = cv_24v96data$humID 
ner_24v96sigIDs = ner_24v96data$humID

x24v96 = list('py'=py_24v96sigIDs,'cv'=cv_24v96sigIDs,'ner'=ner_24v96sigIDs)

venn24v96=Venn(x24v96)
venn24v96

###CHANGE VALUES
fit24v96 <- euler(c('A'=1359,'B'=189,'A&B'=154,'C'=98,'A&C'=63,'B&C'=14,'A&B&C'=19))

Venn24v96 = plot(fit24v96,
                main='24v96 Pairwise Significant Genes',
                fill=viridis(3),
                fill_opacity = .5,
                auto.key=T,
                quantities=T,
                labels=c('Python','Rattlesnake','Watersnake')
)

Venn24v96


