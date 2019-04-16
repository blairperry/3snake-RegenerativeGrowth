
library(pheatmap)
library(viridis)
library(RColorBrewer)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/ipa_results/_new_AnalysisStrategy_02.12.19/UniqueSig_CoreSig_new02.13.19")

cpa.pval <- read.csv('UniqueVsCore_CPA_pval_v3.csv',header = T,row.names = 1)
colnames(cpa.pval) <- c('Py_PRunique','CV_PRunique','Py_Core','CV_Core','Ner_Core')

my_cols <- colorRampPalette(brewer.pal(9,'Blues'))(100)

pheatmap(cpa.pval,
         cluster_cols = F,cluster_rows = F,
         treeheight_row = 0,
         cellheight = 11,cellwidth = 15,
         color = my_cols,border_color = ifelse(cpa.pval > 1.3,'Black','NA'),scale = 'none',
         gaps_col = 2,gaps_row = c(7,26)
         )


### 24v96


cpa.24v96.pval <- read.csv('24v96_IPAresults/24v94_CoreUnique_CPA_pval.csv',header = T,row.names = 1)
colnames(cpa.24v96.pval) <- c('Py_PRunique','CV_PRunique','Py_Core','CV_Core','Ner_Core')

my_cols <- colorRampPalette(brewer.pal(9,'Blues'))(100)

pheatmap(cpa.24v96.pval,
         cluster_cols = F,cluster_rows = F,
         treeheight_row = 0,
         cellheight = 11,cellwidth = 15,
         color = my_cols,border_color = ifelse(cpa.24v96.pval > 1.3,'Black','NA'),scale = 'none',
         gaps_col = 2
)
