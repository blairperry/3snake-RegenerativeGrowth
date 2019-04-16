
#install.packages('ggalluvial')

library(ggalluvial)
library(dplyr)
library(stringr)
library(gridExtra)

setwd('~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/')

######################################################################################################################################################
##
##   PYTHON
##
######################################################################################################################################################

#
# Data input and prep
#
py.fv1 <- read.csv('STAR_Feb2019/pairwise_results/python_0v24_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)
py.1v4 <- read.csv('STAR_Feb2019/pairwise_results/python_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)

py.fv1 <- py.fv1[order(py.fv1$Row.names),]
py.1v4 <- py.1v4[order(py.1v4$Row.names),]


# Determine Up, Down, NotDE for each pairwise comparison
merged.py <- as.data.frame(py.fv1$Row.names)
merged.py$Fv1DPF <- ifelse(py.fv1$log2FoldChange > 0, 'Up','Down')
merged.py$Fv1DPF <- ifelse(py.fv1$IHW_pvalue < 0.05, merged.py$Fv1DPF,'NotDE')
merged.py$Onev4DPF <- ifelse(py.1v4$log2FoldChange > 0, 'Up','Down')
merged.py$Onev4DPF <- ifelse(py.1v4$IHW_pvalue < 0.05, merged.py$Onev4DPF,'NotDE')

#Merge to allow for easier counting
merged.py$joined <- paste(merged.py$Fv1DPF,merged.py$Onev4DPF,sep = '_')

#Count frequency of each combination (i.e. # Up-Up, # Up-Down, etc.)
sigcount.data <- merged.py %>% 
  group_by(joined) %>%
  summarise(Num = length(joined))

#Make alluvial plot input dataframe
py.alluv.data <- as.data.frame(str_split_fixed(sigcount.data$joined, "_", 2))
py.alluv.data$freq <- sigcount.data$Num
py.alluv.data <- py.alluv.data[-c(4,6),]                        # MAY NEED TO UPDATE: remove NotDE-NotDE, NA-NA
py.alluv.data <- py.alluv.data[order(py.alluv.data$V1),]
colnames(py.alluv.data) <- c('Fastedv1DPF','oneDPFv4DPF','Freq')

#
# Alluvial plot
#

py.plot <- ggplot(as.data.frame(py.alluv.data),
       aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/12,reverse=F) +
  geom_stratum(width = 1/12, fill = "grey", color = "black",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Python") +
  ylim(c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())

py.plot


######################################################################################################################################################
##
##   RATTLESNAKE
##
######################################################################################################################################################

#
# Data input and prep
#
cv.fv1 <- read.csv('STAR_Feb2019/pairwise_results/rattlesnake_0v24_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)
cv.1v4 <- read.csv('STAR_Feb2019/pairwise_results/rattlesnake_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)

cv.fv1 <- cv.fv1[order(cv.fv1$Row.names),]
cv.1v4 <- cv.1v4[order(cv.1v4$Row.names),]


# Determine Up, Down, NotDE for each pairwise comparison
merged.cv <- as.data.frame(cv.fv1$Row.names)
merged.cv$Fv1DPF <- ifelse(cv.fv1$log2FoldChange > 0, 'Up','Down')
merged.cv$Fv1DPF <- ifelse(cv.fv1$IHW_pvalue < 0.05, merged.cv$Fv1DPF,'NotDE')
merged.cv$Onev4DPF <- ifelse(cv.1v4$log2FoldChange > 0, 'Up','Down')
merged.cv$Onev4DPF <- ifelse(cv.1v4$IHW_pvalue < 0.05, merged.cv$Onev4DPF,'NotDE')

#Merge to allow for easier counting
merged.cv$joined <- paste(merged.cv$Fv1DPF,merged.cv$Onev4DPF,sep = '_')

#Count frequency of each combination (i.e. # Up-Up, # Up-Down, etc.)
sigcount.data <- merged.cv %>% 
  group_by(joined) %>%
  summarise(Num = length(joined))

#Make alluvial plot input dataframe
cv.alluv.data <- as.data.frame(str_split_fixed(sigcount.data$joined, "_", 2))
cv.alluv.data$freq <- sigcount.data$Num
cv.alluv.data <- cv.alluv.data[-c(4,6),]                        # MAY NEED TO UPDATE: remove NotDE-NotDE, NA-NA
cv.alluv.data <- cv.alluv.data[order(cv.alluv.data$V1),]
colnames(cv.alluv.data) <- c('Fastedv1DPF','oneDPFv4DPF','Freq')

#
# Alluvial plot
#

cv.plot <- ggplot(as.data.frame(cv.alluv.data),
                  aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/12,reverse=F) +
  geom_stratum(width = 1/12, fill = "grey", color = "black",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Rattlesnake") +
  ylim(c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())

cv.plot



######################################################################################################################################################
##
##   WATERSNAKE
##
######################################################################################################################################################


#
# Data input and prep
#
ner.fv1 <- read.csv('STAR_Feb2019/pairwise_results/watersnake_0v24_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)
ner.1v4 <- read.csv('STAR_Feb2019/pairwise_results/watersnake_24v96_DEseq2Out_NonZeroFiltered_02.05.19.csv',header = T)

ner.fv1 <- ner.fv1[order(ner.fv1$Row.names),]
ner.1v4 <- ner.1v4[order(ner.1v4$Row.names),]


# Determine Up, Down, NotDE for each pairwise comparison
merged.ner <- as.data.frame(ner.fv1$Row.names)
merged.ner$Fv1DPF <- ifelse(ner.fv1$log2FoldChange > 0, 'Up','Down')
merged.ner$Fv1DPF <- ifelse(ner.fv1$IHW_pvalue < 0.05, merged.ner$Fv1DPF,'NotDE')
merged.ner$Onev4DPF <- ifelse(ner.1v4$log2FoldChange > 0, 'Up','Down')
merged.ner$Onev4DPF <- ifelse(ner.1v4$IHW_pvalue < 0.05, merged.ner$Onev4DPF,'NotDE')

#Merge to allow for easier counting
merged.ner$joined <- paste(merged.ner$Fv1DPF,merged.ner$Onev4DPF,sep = '_')

#Count frequency of each combination (i.e. # Up-Up, # Up-Down, etc.)
sigcount.data <- merged.ner %>% 
  group_by(joined) %>%
  summarise(Num = length(joined))

#Make alluvial plot input dataframe
ner.alluv.data <- as.data.frame(str_split_fixed(sigcount.data$joined, "_", 2))
ner.alluv.data$freq <- sigcount.data$Num
ner.alluv.data <- ner.alluv.data[-c(3,5),]                        # MAY NEED TO UPDATE: remove NotDE-NotDE, NA-NA
ner.alluv.data <- ner.alluv.data[order(ner.alluv.data$V1),]
colnames(ner.alluv.data) <- c('Fastedv1DPF','oneDPFv4DPF','Freq')

#
# Alluvial plot
#

ner.plot <- ggplot(as.data.frame(ner.alluv.data),
                  aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/12,reverse=F) +
  geom_stratum(width = 1/12, fill = "grey", color = "black",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Watersnake") +
  ylim(c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())

ner.plot

grid.arrange(py.plot,cv.plot,ner.plot,ncol=3)

######################################################################################################################################################
##
##   Format and plot
##
######################################################################################################################################################

py.plot2 <- ggplot(as.data.frame(py.alluv.data),
                  aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/50,reverse=F) +
  geom_stratum(width = 1/50, fill = "grey", color = "white",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  ggtitle("Python") +
  scale_y_continuous(expand = c(0,0),limits=c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())

py.plot2

cv.plot2 <- ggplot(as.data.frame(cv.alluv.data),
                  aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/50,reverse=F) +
  geom_stratum(width = 1/50, fill = "grey", color = "white",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  ggtitle("Rattlesnake") +
  scale_y_continuous(expand = c(0,0),limits=c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())


ner.plot2 <- ggplot(as.data.frame(ner.alluv.data),
                   aes(weight = Freq, axis1 = Fastedv1DPF, axis2 = oneDPFv4DPF)) +
  geom_alluvium(aes(fill = Fastedv1DPF), width = 1/50,reverse=F) +
  geom_stratum(width = 1/50, fill = "grey", color = "white",reverse=F) +
  geom_label(stat = "stratum", label.strata = TRUE,reverse=F) +
  scale_x_continuous(breaks = 1:2, labels = c("Fasted vs. 1DPF", "1DPF vs. 4DPF")) +
  scale_fill_brewer(type = "qual", palette = "Dark2") +
  ggtitle("Watersnake") +
  scale_y_continuous(expand = c(0,0),limits=c(0,3300)) +
  theme(legend.position = 'none',axis.line = element_line(colour = "black"),panel.background = element_blank())


grid.arrange(py.plot2,cv.plot2,ner.plot2,ncol=3)


