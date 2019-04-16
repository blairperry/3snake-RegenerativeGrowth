
library(dplyr)
library(pheatmap)
library(viridis)
library(vegan)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/proteomics/Prot_IPA")

### Main text - 0v1 RNA vs 0v4 Prot

urm.z <- read.csv('./RNAProt_URM_zscore.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.z <- as.data.frame(mutate_all(urm.z, function(x) as.numeric(as.character(x))),row.names=row.names(urm.z))

urm.pval <- read.csv('./RNAProt_URM_pval.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.pval <- as.data.frame(mutate_all(urm.pval, function(x) as.numeric(as.character(x))),row.names=row.names(urm.pval))

urm.z <- urm.z[order(-abs(rowSums(urm.z[,c(1:2)],na.rm=TRUE))),]

urm.pval <- urm.pval[order(match(row.names(urm.pval), row.names(urm.z))), ]

urm.z.top50 <- urm.z[c(0:50),]
urm.pval.top50 <- urm.pval[c(0:50),]

#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(urm.z.top50,na.rm = T)),abs(max(urm.z.top50,na.rm = T)))
myBreaks_0v24_PRW <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

#heatmap
pheatmap(urm.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.pval > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)
pheatmap(urm.z.top50,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.pval.top50 > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)




### Supp - 0v1 and 1v4 RNA vs 0v4 Prot 

urm.full.z <- read.csv('./Full_ProtRNAcomparison_URM_zscore.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.full.z <- as.data.frame(mutate_all(urm.full.z, function(x) as.numeric(as.character(x))),row.names=row.names(urm.full.z))

urm.full.pval <- read.csv('./Full_ProtRNAcomparison_URM_pval.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.full.pval <- as.data.frame(mutate_all(urm.full.pval, function(x) as.numeric(as.character(x))),row.names=row.names(urm.full.pval))

urm.full.z <- urm.full.z[order(-abs(rowSums(urm.full.z[,c(1:2)],na.rm=TRUE))),]

urm.full.pval <- urm.full.pval[order(match(row.names(urm.full.pval), row.names(urm.full.z))), ]

urm.full.z.top50 <- urm.full.z[c(0:50),]
urm.full.pval.top50 <- urm.full.pval[c(0:50),]

#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(urm.full.z.top50,na.rm = T)),abs(max(urm.full.z.top50,na.rm = T)))
myBreaks_0v24_PRW <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

#heatmap
pheatmap(urm.full.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.full.pval > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)
pheatmap(urm.full.z.top50,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.full.pval.top50 > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)



cpa.full.z <- read.csv('./Full_ProtRNAcomparison_cpa_zscore.csv',header=T,row.names=1,stringsAsFactors=FALSE)
cpa.full.z <- as.data.frame(mutate_all(cpa.full.z, function(x) as.numeric(as.character(x))),row.names=row.names(cpa.full.z))

cpa.full.pval <- read.csv('./Full_ProtRNAcomparison_cpa_pval.csv',header=T,row.names=1,stringsAsFactors=FALSE)
cpa.full.pval <- as.data.frame(mutate_all(cpa.full.pval, function(x) as.numeric(as.character(x))),row.names=row.names(cpa.full.pval))

cpa.full.z <- cpa.full.z[order(-abs(rowSums(cpa.full.z[,c(1:2)],na.rm=TRUE))),]

cpa.full.pval <- cpa.full.pval[order(match(row.names(cpa.full.pval), row.names(cpa.full.z))), ]

cpa.full.z.top50 <- cpa.full.z[c(0:50),]
cpa.full.pval.top50 <- cpa.full.pval[c(0:50),]

#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(cpa.full.z.top50,na.rm = T)),abs(max(cpa.full.z.top50,na.rm = T)))
myBreaks_0v24_PRW <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

#heatmap
pheatmap(cpa.full.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(cpa.full.pval > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)
pheatmap(cpa.full.z.top50,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(cpa.full.pval.top50 > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)







#______________________________________________

###
### Split into heatmap based on pattern
###
 
split.z <- cpa.full.z  #change this based on which heatmap needs to be split
split.p <- cpa.full.pval
split.p[is.na(split.p)] <- 0


#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(split.z,na.rm = T)),abs(max(split.z,na.rm = T)))
split_breaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

z.AbsMin = 0   #Set absolute value of minimum z-score

#Pattern1 <- sig. and non-zero all
split.patt1.p <- split.p[which(split.p[,3] > 1.3 & split.p[,1] > 1.3),]
split.patt1.z <- split.z[which(row.names(split.z) %in% row.names(split.patt1.p)),]
#split.patt1.z <- split.patt1.z[which(abs(split.patt1.z[,1]) > 0 & abs(split.patt1.z[,2]) > 0),]
split.patt1.z <- split.patt1.z[which(abs(split.patt1.z[,1]) > z.AbsMin & abs(split.patt1.z[,2]) > z.AbsMin),]
split.patt1.z <- split.patt1.z[order(split.patt1.z[,1] + split.patt1.z[,2],decreasing = T),]
split.patt1.p <- split.patt1.p[which(row.names(split.patt1.p) %in% row.names(split.patt1.z)),]
split.patt1.p <- split.patt1.p[order(match(row.names(split.patt1.p), row.names(split.patt1.z))), ]
colnames(split.patt1.z) <- c('P_RNA','W_RNA','P_Prot','W_Prot')
pheatmap(split.patt1.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt1.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,cellwidth = 15,cellheight = 11,gaps_col = 2,show_colnames = T)



## Supp

#Pattern1 <- sig. in python protein and RNA (and non-zero RNA), not sig in Watersnake protein and RNA
split.patt1.p <- split.p[which(split.p[,1] > 1.3 & split.p[,2] < 1.3 & (split.p[,3] > 1.3 | split.p[,5] > 1.3) & (split.p[,4] < 1.3 & split.p[,6] < 1.3)),]
split.patt1.z <- split.z[which(row.names(split.z) %in% row.names(split.patt1.p)),]
split.patt1.z <- split.patt1.z[which(abs(split.patt1.z[,3]) > 0),]
split.patt1.p <- split.patt1.p[which(row.names(split.patt1.p) %in% row.names(split.patt1.z)),]
split.patt1.p <- split.patt1.p[order(match(row.names(split.patt1.p), row.names(split.patt1.z))), ]
colnames(split.patt1.z) <- c('P.prot','W.prot','P.rnaFv1','W.rnaFv1','P.rna1v4','W.rna1v4')
pheatmap(split.patt1.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt1.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(2,4),show_colnames = T)

#Pattern2 <- sig.in all species and non-zero z-score in at least one
split.patt2.p <- split.p[which(split.p[,1] > 1.3 & split.p[,2] > 1.3 & split.p[,3] > 1.3 & split.p[,4] > 1.3 & split.p[,5] > 1.3 & split.p[,6] > 1.3),]
split.patt2.z <- split.z[which(row.names(split.z) %in% row.names(split.patt2.p)),]
split.patt2.z <- split.patt2.z[which(rowSums(split.patt2.z,na.rm = T) != 0),]
split.patt2.p <- split.patt2.p[which(row.names(split.patt2.p) %in% row.names(split.patt2.z)),]
split.patt2.p <- split.patt2.p[order(match(row.names(split.patt2.p), row.names(split.patt2.z))), ]
colnames(split.patt2.z) <- c('P.prot','W.prot','P.rnaFv1','W.rnaFv1','P.rna1v4','W.rna1v4')
pheatmap(split.patt2.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt2.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(2,4),show_colnames = T)
