
library(dplyr)
library(pheatmap)
library(viridis)
library(vegan)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/ipa_results/_new_AnalysisStrategy_02.12.19/UniqueSig_CoreSig_new02.13.19")


##
## Pairwise Sig Overlap URM heatmap - Fasted vs. 1DPF
##

urm.z <- read.csv('0v24_IPAresults/UniqueVsCore_URM_zscore.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.z <- as.data.frame(mutate_all(urm.z, function(x) as.numeric(as.character(x))),row.names=row.names(urm.z))

urm.pval <- read.csv('0v24_IPAresults/UniqueVsCore_URM_pval.csv',header=T,row.names=1,stringsAsFactors=FALSE)
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




##
## Pairwise Sig Overlap URM heatmap - 1DPF 
##

urm.24v96.z <- read.csv('24v96_IPAresults/24v94_CoreUnique_URM_zscore.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.24v96.z <- as.data.frame(mutate_all(urm.24v96.z, function(x) as.numeric(as.character(x))),row.names=row.names(urm.24v96.z))

urm.24v96.pval <- read.csv('24v96_IPAresults/24v94_CoreUnique_URM_pval.csv',header=T,row.names=1,stringsAsFactors=FALSE)
urm.24v96.pval <- as.data.frame(mutate_all(urm.24v96.pval, function(x) as.numeric(as.character(x))),row.names=row.names(urm.24v96.pval))

urm.24v96.z <- urm.24v96.z[order(-abs(rowSums(urm.24v96.z[,c(1:2)],na.rm=TRUE))),]

urm.24v96.pval <- urm.24v96.pval[order(match(row.names(urm.24v96.pval), row.names(urm.24v96.z))), ]

urm.24v96.z.top50 <- urm.24v96.z[c(0:50),]
urm.24v96.pval.top50 <- urm.24v96.pval[c(0:50),]

#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(urm.24v96.z.top50,na.rm = T)),abs(max(urm.24v96.z.top50,na.rm = T)))
myBreaks_0v24_PRW <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                       seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

#heatmap
pheatmap(urm.24v96.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.24v96.pval > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)
pheatmap(urm.24v96.z.top50,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(urm.24v96.pval.top50 > 1.3,'black',NA),col=mycols,breaks=myBreaks_0v24_PRW,treeheight_row = 0,fontsize_row = 10,cellwidth = 10,cellheight = 10,gaps_col = c(1,2),show_colnames = F)





#______________________________________________

###
### Split into heatmap based on pattern
###
 
split.z <- CPA.z  #change this based on which heatmap needs to be split
split.p <- CPA.pval
split.p[is.na(split.p)] <- 0


#Make symmetrical color pallette
paletteLength <- 40
mycols <- colorRampPalette(c('#136364','#60C2AF','grey95','#EAA155','#D9442E'))(paletteLength)
lim <-  max(abs(min(split.z,na.rm = T)),abs(max(split.z,na.rm = T)))
split_breaks <- c(seq(-(lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                  seq(lim/paletteLength, lim, length.out=floor(paletteLength/2)))

z.AbsMin = 1   #Set absolute value of minimum z-score

#Pattern1 <- sig. and non-zero in PRunique but not Core
split.patt1.p <- split.p[which(split.p[,1] > 1.3 & split.p[,2] > 1.3 & split.p[,3] < 1.3),]
split.patt1.z <- split.z[which(row.names(split.z) %in% row.names(split.patt1.p)),]
#split.patt1.z <- split.patt1.z[which(abs(split.patt1.z[,1]) > 0 & abs(split.patt1.z[,2]) > 0),]
split.patt1.z <- split.patt1.z[which(abs(split.patt1.z[,1]) > z.AbsMin & abs(split.patt1.z[,2]) > z.AbsMin),]
split.patt1.z <- split.patt1.z[order(split.patt1.z[,1] + split.patt1.z[,2],decreasing = T),]
split.patt1.p <- split.patt1.p[which(row.names(split.patt1.p) %in% row.names(split.patt1.z)),]
split.patt1.p <- split.patt1.p[order(match(row.names(split.patt1.p), row.names(split.patt1.z))), ]
colnames(split.patt1.z) <- c('P','R','Pc','Rc','Wc')
pheatmap(split.patt1.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt1.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,cellwidth = 15,cellheight = 11,gaps_col = 2,show_colnames = T)

#Pattern2 <- sig. and non-zero z-score in all species and datasets
split.patt2.p <- split.p[which(split.p[,1] > 1.3 & split.p[,2] > 1.3 & split.p[,3] > 1.3 & split.p[,4] > 1.3 & split.p[,5] > 1.3),]
split.patt2.z <- split.z[which(row.names(split.z) %in% row.names(split.patt2.p)),]
#split.patt2.z <- split.patt2.z[which(abs(split.patt2.z[,1]) != 0 & abs(split.patt2.z[,2]) != 0 & split.patt2.z[,3] != 0),]
split.patt2.z <- split.patt2.z[which(abs(split.patt2.z[,1]) > z.AbsMin & abs(split.patt2.z[,2]) > z.AbsMin & abs(split.patt2.z[,3]) > z.AbsMin & abs(split.patt2.z[,4]) > z.AbsMin & abs(split.patt2.z[,5]) > z.AbsMin),]
split.patt2.z <- split.patt2.z[order(split.patt2.z[,1] + split.patt2.z[,2], decreasing = T),]
split.patt2.p <- split.patt2.p[order(match(row.names(split.patt2.p), row.names(split.patt2.z))), ]
colnames(split.patt2.z) <- c('P','R','Pc','Rc','Wc')
pheatmap(split.patt2.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt2.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,cellwidth = 15,cellheight = 11,gaps_col = 2,show_colnames = T)

#Pattern3 <- sig. and non-zero in Core but not PRunique
split.patt3.p <- split.p[which(split.p[,1] < 1.3 & split.p[,2] < 1.3 & split.p[,3] > 1.3),]
split.patt3.z <- split.z[which(row.names(split.z) %in% row.names(split.patt3.p)),]
#split.patt3.z <- split.patt3.z[which(abs(split.patt3.z[,1]) > 0 & abs(split.patt3.z[,2]) > 0),]
split.patt3.z <- split.patt3.z[which(abs(split.patt3.z[,3]) > z.AbsMin & abs(split.patt3.z[,4]) > z.AbsMin & abs(split.patt3.z[,5]) > z.AbsMin),]
split.patt3.z <- split.patt3.z[order(split.patt3.z[,3] + split.patt3.z[,4],decreasing = T),]
split.patt3.p <- split.patt3.p[which(row.names(split.patt3.p) %in% row.names(split.patt3.z)),]
split.patt3.p <- split.patt3.p[order(match(row.names(split.patt3.p), row.names(split.patt3.z))), ]
colnames(split.patt3.z) <- c('P','R','Pc','Rc','Wc')
pheatmap(split.patt3.z,cluster_cols = F,cluster_rows = F,scale='none',border_color = ifelse(split.patt3.p > 1.3,'black',NA),col=mycols,breaks=split_breaks,treeheight_row = 0,cellwidth = 15,cellheight = 11,gaps_col = 2,show_colnames = T)


 
#______________________________________________________________________________






