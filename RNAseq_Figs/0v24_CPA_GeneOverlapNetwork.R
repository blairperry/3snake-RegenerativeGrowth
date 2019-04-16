#install.packages('GGally')
#install.packages('network')
#install.packages('sna')


library(GGally)
library(network)
library(sna)
library(ggplot2)
library(igraph)
library(viridis)
library(stringr)
library(reshape2)

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/ipa_results/_new_AnalysisStrategy_02.12.19/UniqueSig_CoreSig_new02.13.19/GeneOverlapNetwork")
'%!in%' <- function(x,y)!('%in%'(x,y))

sigPaths <- read.table('0v24_CPA_FilterList.txt',stringsAsFactors = F)

unique <- read.csv('IPA_Output/0v24_PRunique_CPAgeneInfo.csv',header = T,stringsAsFactors = F)
core <- read.csv('IPA_Output/0v24_Core_CPAgeneInfo.csv',header = T,stringsAsFactors = F)


unique.sig <- unique[which(unique$Ingenuity_Canonical_Pathways %in% sigPaths$V1 & unique$X.log.p.value. > 1.3),]
core.sig <- core[which(core$Ingenuity_Canonical_Pathways %in% sigPaths$V1 & core$X.log.p.value. > 1.3),]

shared.sig <- unique.sig[which(unique.sig$Ingenuity_Canonical_Pathways %in% core.sig$Ingenuity_Canonical_Pathways),]
#write.csv(shared.sig[,-6],'Shared_GeneTable_03.12.19.csv')

unique.sig$genes <- strsplit(as.character(unique.sig$Molecules), ",")
unique.sig$NumGenes <- str_count(unique.sig$Molecules,",") + 1
unique.sig$Expanded <- ifelse(unique.sig$Ingenuity_Canonical_Pathways %in% core.sig$Ingenuity_Canonical_Pathways,'Yes','No')

core.sig$genes <- strsplit(as.character(core.sig$Molecules), ",")
core.sig$NumGenes <- str_count(core.sig$Molecules,",") + 1

shared.sig$genes <- strsplit(as.character(shared.sig$Molecules), ",")
shared.sig$NumGenes <- str_count(shared.sig$Molecules,",") + 1

unique.only <- unique.sig[which(unique.sig$Ingenuity_Canonical_Pathways %!in% core.sig$Ingenuity_Canonical_Pathways),]
#write.csv(unique.only[,-6],'UNiqueOnly_GeneTable_03.12.19.csv')

core.only <- core.sig[which(core.sig$Ingenuity_Canonical_Pathways %!in% unique.sig$Ingenuity_Canonical_Pathways),]


### Regen Expanded

cpa.list <- unique.sig$genes

names(cpa.list) <- unique.sig$Ingenuity_Canonical_Pathways


listIntersect <- function(inList) {
  X <- crossprod(table(stack(inList)))
  X[lower.tri(X)] <- NA
  diag(X) <- NA
  out <- na.omit(data.frame(as.table(X)))
  out[order(out$ind), ]
}


overlap.table <- listIntersect(cpa.list)

overlap.table <- merge(overlap.table,unique.sig[,c(1,7)],by.x=1,by.y=1,all.x = T)
overlap.table <- merge(overlap.table,unique.sig[,c(1,7)],by.x=2,by.y=1,all.x = T)

overlap.table$PercentA <- overlap.table$Freq / overlap.table$NumGenes.x *100
overlap.table$PercentB <- overlap.table$Freq / overlap.table$NumGenes.y *100

overlap.table$MaxPerc <- round(pmax(overlap.table$PercentA,overlap.table$PercentB)) / 100


min.overlap.thresh <- 0.5 #minimum number of connections to warrant an edge in network

filtered.overlap.table <- overlap.table[which(overlap.table$MaxPerc >= min.overlap.thresh),]



g <- graph.data.frame(filtered.overlap.table, directed=F)

net <- as.matrix(get.adjacency(g,attr='MaxPerc'))
net <- network(net,directed = F,ignore.eval=F,names.eval='MaxPerc')



cpa.net <- ggnet2(net,
                  shape=16,
                  palette = 'Set2',
                  mode= 'kamadakawai', layout.par = list(niter = 10000), layout.exp = .5,
                  label=T, label.size = 3,
                  #edge.label = 'MaxPerc', edge.label.size = 2,
                  edge.size = 'MaxPerc', edge.alpha = .5,
                  size = "indegree"
)

cpa.net



### Regen Unique

cpa.list <- unique.only$genes

names(cpa.list) <- unique.only$Ingenuity_Canonical_Pathways


listIntersect <- function(inList) {
  X <- crossprod(table(stack(inList)))
  X[lower.tri(X)] <- NA
  diag(X) <- NA
  out <- na.omit(data.frame(as.table(X)))
  out[order(out$ind), ]
}


overlap.table <- listIntersect(cpa.list)

overlap.table <- merge(overlap.table,unique.sig[,c(1,7)],by.x=1,by.y=1,all.x = T)
overlap.table <- merge(overlap.table,unique.sig[,c(1,7)],by.x=2,by.y=1,all.x = T)

overlap.table$PercentA <- overlap.table$Freq / overlap.table$NumGenes.x *100
overlap.table$PercentB <- overlap.table$Freq / overlap.table$NumGenes.y *100

overlap.table$MaxPerc <- round(pmax(overlap.table$PercentA,overlap.table$PercentB)) / 100


min.overlap.thresh <- 0.5 #minimum number of connections to warrant an edge in network

filtered.overlap.table <- overlap.table[which(overlap.table$MaxPerc >= min.overlap.thresh),]



g <- graph.data.frame(filtered.overlap.table, directed=F)

net <- as.matrix(get.adjacency(g,attr='MaxPerc'))
net <- network(net,directed = F,ignore.eval=F,names.eval='MaxPerc')


cpa.net <- ggnet2(net,
                  shape=16,
                  palette = 'Set2',
                  mode= 'kamadakawai', layout.par = list(niter = 10000), layout.exp = .5,
                  label=T, label.size = 3,
                  #edge.label = 'MaxPerc', edge.label.size = 2,
                  edge.size = 'MaxPerc', edge.alpha = .5,
                  size = "indegree"
)

cpa.net



### Core

cpa.list <- core.sig$genes

names(cpa.list) <- core.sig$Ingenuity_Canonical_Pathways


listIntersect <- function(inList) {
  X <- crossprod(table(stack(inList)))
  X[lower.tri(X)] <- NA
  diag(X) <- NA
  out <- na.omit(data.frame(as.table(X)))
  out[order(out$ind), ]
}


overlap.table <- listIntersect(cpa.list)

overlap.table <- merge(overlap.table,core.sig[,c(1,7)],by.x=1,by.y=1,all.x = T)
overlap.table <- merge(overlap.table,core.sig[,c(1,7)],by.x=2,by.y=1,all.x = T)

overlap.table$PercentA <- overlap.table$Freq / overlap.table$NumGenes.x *100
overlap.table$PercentB <- overlap.table$Freq / overlap.table$NumGenes.y *100

overlap.table$MaxPerc <- round(pmax(overlap.table$PercentA,overlap.table$PercentB)) / 100


min.overlap.thresh <- 0.5 #minimum number of connections to warrant an edge in network

filtered.overlap.table <- overlap.table[which(overlap.table$MaxPerc >= min.overlap.thresh),]



g <- graph.data.frame(filtered.overlap.table, directed=F)

net <- as.matrix(get.adjacency(g,attr='MaxPerc'))
net <- network(net,directed = F,ignore.eval=F,names.eval='MaxPerc')


cpa.net <- ggnet2(net,
                  shape=16,
                  palette = 'Set2',
                  mode= 'kamadakawai', layout.par = list(niter = 10000), layout.exp = .5,
                  label=T, label.size = 3,
                  #edge.label = 'MaxPerc', edge.label.size = 2,
                  edge.size = 'MaxPerc', edge.alpha = .5,
                  size = "indegree"
)

cpa.net

