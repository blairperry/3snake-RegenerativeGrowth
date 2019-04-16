
###
### RNA Abundance vs. Protein Abundance 
###

setwd("~/Dropbox/CastoeLabFolder/projects/SnakePhysiolRemod/__3Species_Intestine_MS/data/proteomics")

# Get average normPSM for each time point

Py.normPSM <- read.csv('Python_IDConvert/Py.DEseq2NormPSM_11Apr18_IDconvert_03.22.19.csv',stringsAsFactors = F,header=T)

Py.AvgPSM <- data.frame(Py.normPSM$X.3)
Py.AvgPSM$AvgPSM.0 <- rowMeans(Py.normPSM[,c(5,6,7)])
Py.AvgPSM$AvgPSM.24 <- rowMeans(Py.normPSM[,c(8,9,10)])
Py.AvgPSM$AvgPSM.96 <- rowMeans(Py.normPSM[,c(11,12)])
Py.AvgPSM <- Py.AvgPSM[which(!is.na(Py.AvgPSM$Py.normPSM.X.3)),]

## Read in norm count file in order to calc. average norm counts per timepoint

Py.RNA.NormCounts <- read.csv('../STAR_Feb2019/norm_counts/py_SI_normcounts_DEseq_02.05.2019.csv',stringsAsFactors = F,header=T)

Py.RNA.AvgNormCount <- data.frame(Py.RNA.NormCounts$X)
Py.RNA.AvgNormCount$AvgNormCount.0 <- rowMeans(Py.RNA.NormCounts[,c(2,3,4,5,6)])
Py.RNA.AvgNormCount$AvgNormCount.24 <- rowMeans(Py.RNA.NormCounts[,c(7,8,9,10,11,12)])
Py.RNA.AvgNormCount$AvgNormCount.96 <- rowMeans(Py.RNA.NormCounts[,c(13,14,15,16,17)])

pyNum_to_hum <- read.table('Python_IDConvert/PyGeneNum_to_Human.txt')

Py.RNA.AvgNormCount.humID <- merge(Py.RNA.AvgNormCount,pyNum_to_hum,by.x=1,by.y=1,all.x=T)
Py.RNA.AvgNormCount.humID <- Py.RNA.AvgNormCount.humID[which(!is.na(Py.RNA.AvgNormCount.humID$V2)),]

## Fasted 

## Merge protein and RNA data

Py.AvgPSM.0 <- Py.AvgPSM[,c(1,2)]
Py.AvgPSM.0 <- Py.AvgPSM.0[which(Py.AvgPSM.0$AvgPSM.0 > 0),]

Py.RNA.AvgNormCount.0 <- Py.RNA.AvgNormCount.humID[,c(5,2)]


Py.RNAProt.AvgAbund.0 <- merge(Py.AvgPSM.0,Py.RNA.AvgNormCount.0,by.x=1,by.y = 1,all=F)

plot(Py.RNAProt.AvgAbund.0$AvgNormCount.0,Py.RNAProt.AvgAbund.0$AvgPSM.0,pch=19,col=rgb(0,0.4,0.4,alpha=0.4),main='RNA - Fasted vs 1DPF\n Protein - Fasted vs 1DPF')

## Log2 Transform

Py.RNAProt.AvgAbund.0.log2 <- Py.RNAProt.AvgAbund.0
Py.RNAProt.AvgAbund.0.log2$AvgPSM.0 <- log2(Py.RNAProt.AvgAbund.0.log2$AvgPSM.0)
Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0 <- log2(Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0)
Py.RNAProt.AvgAbund.0.log2 <- Py.RNAProt.AvgAbund.0.log2[which(!is.infinite(Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0)),]

par(pty='s')
plot(Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0,Py.RNAProt.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - Fasted',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(-5,15),ylim=c(-2,7.5))

# Linear regression on log transformed data
AvgAbund.0.lm <- lm(formula = AvgPSM.0 ~ AvgNormCount.0, data = Py.RNAProt.AvgAbund.0.log2)
summary(AvgAbund.0.lm)
abline(b=0.10541,a=1.52247)
text(x=0,y=7,labels = expression(R^2 == 0.02965))
text(x=0,y=6.5,labels='p < 0.001')

#Subset RNA min 5 and Prot min 0
Py.RNAmin5.Protmin0.AvgAbund.0.log2 <- Py.RNAProt.AvgAbund.0.log2[which(Py.RNAProt.AvgAbund.0.log2$AvgPSM.0 > 0 &  Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0 > 5),]

Py.RNAmin5.Protmin0.AvgAbund.0.log2.exclude <- Py.RNAProt.AvgAbund.0.log2[!(Py.RNAProt.AvgAbund.0.log2$Py.normPSM.X %in% Py.RNAmin5.Protmin0.AvgAbund.0.log2$Py.normPSM.X),]

plot(Py.RNAmin5.Protmin0.AvgAbund.0.log2.exclude$AvgNormCount.0,Py.RNAmin5.Protmin0.AvgAbund.0.log2.exclude$AvgPSM.0,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - Fasted',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,15),ylim=c(-2,7.5))
points(Py.RNAmin5.Protmin0.AvgAbund.0.log2$AvgNormCount.0,Py.RNAmin5.Protmin0.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(0.21,0.28,0.51,alpha=0.6))

#all data regression line
abline(b=0.10541,a=1.52247,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.02965),col='grey40')
text(x=1,y=5,labels='p < 0.001',col='grey40')

RNAmin5.Protmin0.AvgAbund.0.lm <- lm(formula = AvgPSM.0 ~ AvgNormCount.0, data = Py.RNAmin5.Protmin0.AvgAbund.0.log2)
summary(RNAmin5.Protmin0.AvgAbund.0.lm)
abline(b=0.17189,a=1.12714)
text(x=1,y=7,labels = expression(R^2 == 0.05143))
text(x=1,y=6.5,labels='p =2.178e-05')


## 1DPF

## Merge protein and RNA data

Py.AvgPSM.24 <- Py.AvgPSM[,c(1,3)]
Py.AvgPSM.24 <- Py.AvgPSM.24[which(Py.AvgPSM.24$AvgPSM.24 > 0),]

Py.RNA.AvgNormCount.24 <- Py.RNA.AvgNormCount.humID[,c(5,3)]

Py.RNAProt.AvgAbund.24 <- merge(Py.AvgPSM.24,Py.RNA.AvgNormCount.24,by.x=1,by.y =1,all=F)

#plot(Py.RNAProt.AvgAbund.24$AvgNormCount.24,Py.RNAProt.AvgAbund.24$AvgPSM.24,pch=19,col=rgb(0,0.4,0.4,alpha=0.4),main='RNA - Fasted vs 1DPF\n Protein - Fasted vs 1DPF')

Py.RNAProt.AvgAbund.24.log2 <- Py.RNAProt.AvgAbund.24
Py.RNAProt.AvgAbund.24.log2$AvgPSM.24 <- log2(Py.RNAProt.AvgAbund.24.log2$AvgPSM.24)
Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24 <- log2(Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24)
Py.RNAProt.AvgAbund.24.log2 <- Py.RNAProt.AvgAbund.24.log2[which(!is.infinite(Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24)),]

plot(Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24,Py.RNAProt.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - 1DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(-5,15),ylim=c(-2,7.5))

AvgAbund.24.lm <- lm(formula = AvgPSM.24 ~ AvgNormCount.24, data = Py.RNAProt.AvgAbund.24.log2)
summary(AvgAbund.24.lm)
abline(b=0.14301,a=0.96884)
text(x=0,y=7,labels = expression(R^2 == 0.05149))
text(x=0,y=6.5,labels='p = 1.609e-06')


#Subset RNA min 5 and Prot min 0
Py.RNAmin5.Protmin0.AvgAbund.24.log2 <- Py.RNAProt.AvgAbund.24.log2[which(Py.RNAProt.AvgAbund.24.log2$AvgPSM.24 > 0 &  Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24 > 5),]

Py.RNAmin5.Protmin0.AvgAbund.24.log2.exclude <- Py.RNAProt.AvgAbund.24.log2[!(Py.RNAProt.AvgAbund.24.log2$Py.normPSM.X %in% Py.RNAmin5.Protmin0.AvgAbund.24.log2$Py.normPSM.X),]

plot(Py.RNAmin5.Protmin0.AvgAbund.24.log2.exclude$AvgNormCount.24,Py.RNAmin5.Protmin0.AvgAbund.24.log2.exclude$AvgPSM.24,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - 1DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,15),ylim=c(-2,7.5))
points(Py.RNAmin5.Protmin0.AvgAbund.24.log2$AvgNormCount.24,Py.RNAmin5.Protmin0.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(0.21,0.28,0.51,alpha=0.6))

#all data regression line
abline(b=0.14301,a=0.96884,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.05149),col='grey40')
text(x=1,y=5,labels='p = 1.609e-06',col='grey40')

RNAmin5.Protmin0.AvgAbund.24.lm <- lm(formula = AvgPSM.24 ~ AvgNormCount.24, data = Py.RNAmin5.Protmin0.AvgAbund.24.log2)
summary(RNAmin5.Protmin0.AvgAbund.24.lm)
abline(b=0.42091,a=-1.02834)
text(x=1,y=7,labels = expression(R^2 == 0.2013))
text(x=1,y=6.5,labels='p < 2.2e-16')




## 4DPF

## Merge protein and RNA data

Py.AvgPSM.96 <- Py.AvgPSM[,c(1,4)]
Py.AvgPSM.96 <- Py.AvgPSM.96[which(Py.AvgPSM.96$AvgPSM.96 > 0),]

Py.RNA.AvgNormCount.96 <- Py.RNA.AvgNormCount.humID[,c(5,4)]

Py.RNAProt.AvgAbund.96 <- merge(Py.AvgPSM.96,Py.RNA.AvgNormCount.96,by.x=1,by.y = 1,all=F)

#plot(Py.RNAProt.AvgAbund.96$AvgNormCount.96,Py.RNAProt.AvgAbund.96$AvgPSM.96,pch=19,col=rgb(0,0.4,0.4,alpha=0.4),main='RNA - Fasted vs 1DPF\n Protein - Fasted vs 1DPF')

Py.RNAProt.AvgAbund.96.log2 <- Py.RNAProt.AvgAbund.96
Py.RNAProt.AvgAbund.96.log2$AvgPSM.96 <- log2(Py.RNAProt.AvgAbund.96.log2$AvgPSM.96)
Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96 <- log2(Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96)
Py.RNAProt.AvgAbund.96.log2 <- Py.RNAProt.AvgAbund.96.log2[which(!is.infinite(Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96)),]

plot(Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96,Py.RNAProt.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - 4DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(-5,15),ylim=c(-2,7.5))

AvgAbund.96.lm <- lm(formula = AvgPSM.96 ~ AvgNormCount.96, data = Py.RNAProt.AvgAbund.96.log2)
summary(AvgAbund.96.lm)
abline(b=0.10077,a=1.49432)
text(x=0,y=7,labels = expression(R^2 == 0.01433))
text(x=0,y=6.5,labels='p = 0.01592')



#Subset RNA min 5 and Prot min 0
Py.RNAmin5.Protmin0.AvgAbund.96.log2 <- Py.RNAProt.AvgAbund.96.log2[which(Py.RNAProt.AvgAbund.96.log2$AvgPSM.96 > 0 &  Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96 > 5),]

Py.RNAmin5.Protmin0.AvgAbund.96.log2.exclude <- Py.RNAProt.AvgAbund.96.log2[!(Py.RNAProt.AvgAbund.96.log2$Py.normPSM.X %in% Py.RNAmin5.Protmin0.AvgAbund.96.log2$Py.normPSM.X),]

plot(Py.RNAmin5.Protmin0.AvgAbund.96.log2.exclude$AvgNormCount.96,Py.RNAmin5.Protmin0.AvgAbund.96.log2.exclude$AvgPSM.96,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - 4DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,15),ylim=c(-2,7.5))
points(Py.RNAmin5.Protmin0.AvgAbund.96.log2$AvgNormCount.96,Py.RNAmin5.Protmin0.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(0.21,0.28,0.51,alpha=0.6))

#all data regression line
abline(b=0.10077,a=1.49432,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.01433),col='grey40')
text(x=1,y=5,labels='p = 0.01592',col='grey40')

RNAmin5.Protmin0.AvgAbund.96.lm <- lm(formula = AvgPSM.96 ~ AvgNormCount.96, data = Py.RNAmin5.Protmin0.AvgAbund.96.log2)
summary(RNAmin5.Protmin0.AvgAbund.96.lm)
abline(b=0.42479,a=-0.83445)
text(x=1,y=7,labels = expression(R^2 == 0.1455))
text(x=1,y=6.5,labels='p = 9.281e-11')



par(mfrow=c(3,1))

#plot all three

par(mfrow=c(1,1))


###
###
### WATERSNAKE
###
###


###
### RNA Abundance vs. Protein Abundance
###

# Get average normPSM for each time point

Ner.normPSM <- read.csv('PSM/Ner.DEseq2NormPSM_11Apr18.csv',header=T)

Ner.AvgPSM <- data.frame(Ner.normPSM$Ner.Tsirt._ID)
Ner.AvgPSM$AvgPSM.0 <- rowMeans(Ner.normPSM[,c(3,4,5)])
Ner.AvgPSM$AvgPSM.24 <- rowMeans(Ner.normPSM[,c(6,7,8)])
Ner.AvgPSM$AvgPSM.96 <- rowMeans(Ner.normPSM[,c(9,10,11)])


## Read in norm count file in order to calc. average norm counts per timepoint

Ner.RNA.NormCounts <- read.csv('../STAR_Feb2019/norm_counts/ner_SI_normcounts_DEseq_02.05.2019.csv',stringsAsFactors = F,header=T)

Ner.RNA.AvgNormCount <- data.frame(Ner.RNA.NormCounts$X)
Ner.RNA.AvgNormCount$AvgNormCount.0 <- rowMeans(Ner.RNA.NormCounts[,c(2,3,4,5)])
Ner.RNA.AvgNormCount$AvgNormCount.24 <- rowMeans(Ner.RNA.NormCounts[,c(6,7,8,9)])
Ner.RNA.AvgNormCount$AvgNormCount.96 <- rowMeans(Ner.RNA.NormCounts[,c(10,11,12,13)])

NerNum_to_Tsirt <- read.table('Python_IDConvert/NerGeneNum_to_Tsirt.txt')

Ner.RNA.AvgNormCount.HumID <- merge(Ner.RNA.AvgNormCount,NerNum_to_Tsirt,by.x=1,by.y=1,all.x=T)
Ner.RNA.AvgNormCount.HumID <- Ner.RNA.AvgNormCount.HumID[which(!is.na(Ner.RNA.AvgNormCount.HumID$V2)),]

## Fasted

## Merge protein and RNA data

Ner.AvgPSM.0 <- Ner.AvgPSM[,c(1,2)]
Ner.AvgPSM.0 <- Ner.AvgPSM.0[which(Ner.AvgPSM.0$AvgPSM.0 > 0),]

Ner.RNA.AvgNormCount.0 <- Ner.RNA.AvgNormCount.HumID[,c(5,2)]

Ner.RNAProt.AvgAbund.0 <- merge(Ner.AvgPSM.0,Ner.RNA.AvgNormCount.0,by.x=1,by.y = 1,all=F)

## Log2 Transform

Ner.RNAProt.AvgAbund.0.log2 <- Ner.RNAProt.AvgAbund.0
Ner.RNAProt.AvgAbund.0.log2$AvgPSM.0 <- log2(Ner.RNAProt.AvgAbund.0.log2$AvgPSM.0)
Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0 <- log2(Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0)
Ner.RNAProt.AvgAbund.0.log2 <- Ner.RNAProt.AvgAbund.0.log2[which(!is.infinite(Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0)),]

par(pty='s')
plot(Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0,Ner.RNAProt.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - Fasted',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(-2,17),ylim=c(-2,7))

# Linear regression on log transformed data
AvgAbund.0.lm <- lm(formula = AvgPSM.0 ~ AvgNormCount.0, data = Ner.RNAProt.AvgAbund.0.log2)
summary(AvgAbund.0.lm)
abline(b=0.15082,a=0.93832)
text(x=0,y=6,labels = expression(R^2 == 0.09651))
text(x=0,y=5.5,labels='p = 4.202e-14')


#Subset RNA min 5 and Prot min 0
Ner.RNAmin5.Protmin0.AvgAbund.0.log2 <- Ner.RNAProt.AvgAbund.0.log2[which(Ner.RNAProt.AvgAbund.0.log2$AvgPSM.0 > 0 &  Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0 > 5),]

Ner.RNAmin5.Protmin0.AvgAbund.0.log2.exclude <- Ner.RNAProt.AvgAbund.0.log2[!(Ner.RNAProt.AvgAbund.0.log2[,1] %in% Ner.RNAmin5.Protmin0.AvgAbund.0.log2[,1]),]

plot(Ner.RNAmin5.Protmin0.AvgAbund.0.log2.exclude$AvgNormCount.0,Ner.RNAmin5.Protmin0.AvgAbund.0.log2.exclude$AvgPSM.0,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - Fasted',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(-2,20),ylim=c(-2,7.5))
points(Ner.RNAmin5.Protmin0.AvgAbund.0.log2$AvgNormCount.0,Ner.RNAmin5.Protmin0.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(0.55,0.23,0.22,alpha=0.5))

#all data regression line
abline(b=0.15082,a=0.93832,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.09651),col='grey40')
text(x=1,y=5,labels='p = 4.202e-14',col='grey40')

RNAmin5.Protmin0.AvgAbund.0.lm <- lm(formula = AvgPSM.0 ~ AvgNormCount.0, data = Ner.RNAmin5.Protmin0.AvgAbund.0.log2)
summary(RNAmin5.Protmin0.AvgAbund.0.lm)
abline(b=0.20918,a=0.45045)
text(x=1,y=7,labels = expression(R^2 == 0.1266))
text(x=1,y=6.5,labels='p = 8.86e-16')





## 1DPF

## Merge protein and RNA data

Ner.AvgPSM.24 <- Ner.AvgPSM[,c(1,3)]
Ner.AvgPSM.24 <- Ner.AvgPSM.24[which(Ner.AvgPSM.24$AvgPSM.24 > 0),]

Ner.RNA.AvgNormCount.24 <- Ner.RNA.AvgNormCount.HumID[,c(5,3)]

Ner.RNAProt.AvgAbund.24 <- merge(Ner.AvgPSM.24,Ner.RNA.AvgNormCount.24,by.x=1,by.y = 1,all=F)

#plot(Ner.RNAProt.AvgAbund.24$AvgNormCount.24,Ner.RNAProt.AvgAbund.24$AvgPSM.24,pch=19,col=rgb(0,0.4,0.4,alpha=0.4),main='RNA - Fasted vs 1DPF\n Protein - Fasted vs 1DPF')

Ner.RNAProt.AvgAbund.24.log2 <- Ner.RNAProt.AvgAbund.24
Ner.RNAProt.AvgAbund.24.log2$AvgPSM.24 <- log2(Ner.RNAProt.AvgAbund.24.log2$AvgPSM.24)
Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24 <- log2(Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24)
Ner.RNAProt.AvgAbund.24.log2 <- Ner.RNAProt.AvgAbund.24.log2[which(!is.infinite(Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24)),]

plot(Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24,Ner.RNAProt.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - 1DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,20),ylim=c(-2,7.5))


############ STOPPED HERE

AvgAbund.24.lm <- lm(formula = AvgPSM.24 ~ AvgNormCount.24, data = Ner.RNAProt.AvgAbund.24.log2)
summary(AvgAbund.24.lm)
abline(b=0.21060,a=0.06689)
text(x=1,y=7,labels = expression(R^2 == 0.1491))
text(x=1,y=6.5,labels='p < 2.2e-16')


#Subset RNA min 5 and Prot min 0
Ner.RNAmin5.Protmin0.AvgAbund.24.log2 <- Ner.RNAProt.AvgAbund.24.log2[which(Ner.RNAProt.AvgAbund.24.log2$AvgPSM.24 > 0 &  Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24 > 5),]

Ner.RNAmin5.Protmin0.AvgAbund.24.log2.exclude <- Ner.RNAProt.AvgAbund.24.log2[!(Ner.RNAProt.AvgAbund.24.log2[,1] %in% Ner.RNAmin5.Protmin0.AvgAbund.24.log2[,1]),]

plot(Ner.RNAmin5.Protmin0.AvgAbund.24.log2.exclude$AvgNormCount.24,Ner.RNAmin5.Protmin0.AvgAbund.24.log2.exclude$AvgPSM.24,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - 1DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,20),ylim=c(-2,7.5))
points(Ner.RNAmin5.Protmin0.AvgAbund.24.log2$AvgNormCount.24,Ner.RNAmin5.Protmin0.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(0.55,0.23,0.22,alpha=0.5))

#all data regression line
abline(b=0.21060,a=0.06689,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.1491),col='grey40')
text(x=1,y=5,labels='p < 2.2e-16',col='grey40')

RNAmin5.Protmin0.AvgAbund.24.lm <- lm(formula = AvgPSM.24 ~ AvgNormCount.24, data = Ner.RNAmin5.Protmin0.AvgAbund.24.log2)
summary(RNAmin5.Protmin0.AvgAbund.24.lm)
abline(b=0.27713,a=-0.38015)
text(x=1,y=7,labels = expression(R^2 == 0.212))
text(x=1,y=6.5,labels='p < 2.2e-16')




## 4DPF

## Merge protein and RNA data

Ner.AvgPSM.96 <- Ner.AvgPSM[,c(1,4)]
Ner.AvgPSM.96 <- Ner.AvgPSM.96[which(Ner.AvgPSM.96$AvgPSM.96 > 0),]

Ner.RNA.AvgNormCount.96 <- Ner.RNA.AvgNormCount.HumID[,c(5,4)]

Ner.RNAProt.AvgAbund.96 <- merge(Ner.AvgPSM.96,Ner.RNA.AvgNormCount.96,by.x=1,by.y = 1,all=F)

#plot(Ner.RNAProt.AvgAbund.96$AvgNormCount.96,Ner.RNAProt.AvgAbund.96$AvgPSM.96,pch=19,col=rgb(0,0.4,0.4,alpha=0.4),main='RNA - Fasted vs 1DPF\n Protein - Fasted vs 1DPF')

Ner.RNAProt.AvgAbund.96.log2 <- Ner.RNAProt.AvgAbund.96
Ner.RNAProt.AvgAbund.96.log2$AvgPSM.96 <- log2(Ner.RNAProt.AvgAbund.96.log2$AvgPSM.96)
Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96 <- log2(Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96)
Ner.RNAProt.AvgAbund.96.log2 <- Ner.RNAProt.AvgAbund.96.log2[which(!is.infinite(Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96)),]

plot(Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96,Ner.RNAProt.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(0.3,0.6,0.3,alpha=0.5),main='RNA vs. Protein Abundance - 4DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,20),ylim=c(-2,7.5))

AvgAbund.96.lm <- lm(formula = AvgPSM.96 ~ AvgNormCount.96, data = Ner.RNAProt.AvgAbund.96.log2)
summary(AvgAbund.96.lm)
abline(b=0.25767,a=-0.37724)
text(x=0,y=7,labels = expression(R^2 == 0.1174))
text(x=0,y=6.5,labels='p = 2.758e-14')



#Subset RNA min 5 and Prot min 0
Ner.RNAmin5.Protmin0.AvgAbund.96.log2 <- Ner.RNAProt.AvgAbund.96.log2[which(Ner.RNAProt.AvgAbund.96.log2$AvgPSM.96 > 0 &  Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96 > 5),]

Ner.RNAmin5.Protmin0.AvgAbund.96.log2.exclude <- Ner.RNAProt.AvgAbund.96.log2[!(Ner.RNAProt.AvgAbund.96.log2$Ner.normPSM.Ner.Tsirt._ID %in% Ner.RNAmin5.Protmin0.AvgAbund.96.log2$Ner.normPSM.Ner.Tsirt._ID),]

plot(Ner.RNAmin5.Protmin0.AvgAbund.96.log2.exclude$AvgNormCount.96,Ner.RNAmin5.Protmin0.AvgAbund.96.log2.exclude$AvgPSM.96,pch=19,col=rgb(0.6,0.6,0.6,alpha=0.5),main='RNA vs. Protein Abundance - 4DPF',ylab= 'Protein Abundance - Log2(Avg. %PSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)',xlim=c(0,20),ylim=c(-2,7.5))
points(Ner.RNAmin5.Protmin0.AvgAbund.96.log2$AvgNormCount.96,Ner.RNAmin5.Protmin0.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(0.55,0.23,0.22,alpha=0.5))

#all data regression line
abline(b=0.25767,a=-0.37724,col='grey40',lty=2)
text(x=1,y=5.5,labels = expression(R^2 == 0.1174),col='grey40')
text(x=1,y=5,labels='p < 0.001',col='grey40')

RNAmin5.Protmin0.AvgAbund.96.lm <- lm(formula = AvgPSM.96 ~ AvgNormCount.96, data = Ner.RNAmin5.Protmin0.AvgAbund.96.log2)
summary(RNAmin5.Protmin0.AvgAbund.96.lm)
abline(b=0.30232,a=-0.40650)
text(x=1,y=7,labels = expression(R^2 == 0.1706))
text(x=1,y=6.5,labels='p < 2.2e-16')



par(mfrow=c(3,1))

#plot all three

par(mfrow=c(1,1))





# ###
# ### Combine all 3 species into 3x3 figure
# ###
# 
# par(mfrow=c(3,3))
# #python F
# plot(Py.RNAProt.AvgAbund.0.log2$AvgNormCount.0,Py.RNAProt.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(0.38,0.38,0.38,alpha=0.4),main='Python',xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.09900,a=1.51344)
# text(x=0,y=7,labels = expression(R^2 == 0.03029))
# text(x=0,y=6,labels='p < 0.001')
# 
# #rattlesnake F
# plot(CV.RNAProt.AvgAbund.0.log2$AvgNormCount.0,CV.RNAProt.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(.29,.73,.45,alpha=0.4),main='Rattlesnake',xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.02981,a=1.74275)
# text(x=0,y=6,labels = expression(R^2 == 0.002459))
# text(x=0,y=5.5,labels='p = 0.1758')
# 
# #watersnake F
# plot(Ner.RNAProt.AvgAbund.0.log2$AvgNormCount.0,Ner.RNAProt.AvgAbund.0.log2$AvgPSM.0,pch=19,col=rgb(.76,.54,.09,alpha=0.4),main='Watersnake',xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.02981,a=1.74275)
# text(x=1,y=6,labels = expression(R^2 == 0.002459))
# text(x=1,y=5.5,labels='p = 0.1758')
# 
# #python 1
# plot(Py.RNAProt.AvgAbund.24.log2$AvgNormCount.24,Py.RNAProt.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(0.38,0.38,0.38,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.15706,a=0.73743)
# text(x=0,y=7,labels = expression(R^2 == 0.06726))
# text(x=0,y=6,labels='p < 0.001')
# 
# #rattlesnake 1
# plot(CV.RNAProt.AvgAbund.24.log2$AvgNormCount.24,CV.RNAProt.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(.29,.73,.45,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.05322,a=1.18433)
# text(x=0,y=6,labels = expression(R^2 == 0.006782))
# text(x=0,y=5.5,labels='p = 0.08608')
# 
# #watersnake 1
# plot(Ner.RNAProt.AvgAbund.24.log2$AvgNormCount.24,Ner.RNAProt.AvgAbund.24.log2$AvgPSM.24,pch=19,col=rgb(.76,.54,.09,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.13773,a=0.85845)
# text(x=1,y=6,labels = expression(R^2 == 0.07112))
# text(x=1,y=5.5,labels='p < 0.001')
# 
# #python 4
# plot(Py.RNAProt.AvgAbund.96.log2$AvgNormCount.96,Py.RNAProt.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(0.38,0.38,0.38,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.16670,a=0.89903)
# text(x=0,y=7,labels = expression(R^2 == 0.05387))
# text(x=0,y=6,labels='p < 0.001')
# 
# #rattlesnake 4
# plot(CV.RNAProt.AvgAbund.96.log2$AvgNormCount.96,CV.RNAProt.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(.29,.73,.45,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.02632,a=1.14442)
# text(x=0,y=6,labels = expression(R^2 == 0))
# text(x=0,y=5.5,labels='p = 0.4166')
# 
# #watersnake 4
# plot(Ner.RNAProt.AvgAbund.96.log2$AvgNormCount.96,Ner.RNAProt.AvgAbund.96.log2$AvgPSM.96,pch=19,col=rgb(.76,.54,.09,alpha=0.4),xlim=c(-5,17),ylim=c(-2,7.5),ylab= 'Protein Abundance - Log2(Avg. NormPSM)',xlab='RNA Abundance - Log2(Avg. Norm. Counts)')
# abline(b=0.08427,a=1.59341)
# text(x=1,y=6,labels = expression(R^2 == 0.039))
# text(x=1,y=5.5,labels='p < 0.001')
# par(mfrow=c(1,1))
# 
# 
# ###
# ### Combine Fold Change comparisons figures
# ###
# 
# par(mfrow=c(2,3))
# #python a
# plot(Py.Merged.0v24.FC.sig$RNA_Log2FC,Py.Merged.0v24.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 1DPF - Log2FC',col='red',main='Python',xlim=c(-6,6),ylim=c(-6,6))
# points(Py.Merged.0v24.FC.notsig$RNA_Log2FC,Py.Merged.0v24.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.05428,a=-0.06859)
# text(x=-5,y=4.5,labels = expression(R^2 == 0.001329))
# text(x=-5,y=4,labels='p = 0.2254')
# abline(b=0.10952,a=-0.29071,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.008342),col='blue')
# text(x=-5,y=3,labels='p < 0.01',col='blue')
# 
# 
# #rattlesnake a
# plot(CV.Merged.0v24.FC.sig$RNA_Log2FC,CV.Merged.0v24.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 1DPF - Log2FC',col='red',main='Rattlesnake',xlim=c(-6,6),ylim=c(-6,6))
# points(CV.Merged.0v24.FC.notsig$RNA_Log2FC,CV.Merged.0v24.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.25008,a=-0.70342)
# text(x=-5,y=4.5,labels = expression(R^2 == 0.04544))
# text(x=-5,y=4,labels='p < 0.05')
# abline(b=0.25989,a=-0.91478,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.02104),col='blue')
# text(x=-5,y=3,labels='p < 0.01',col='blue')
# 
# 
# #watersnake a
# plot(Ner.Merged.0v24.FC.sig$RNA_Log2FC,Ner.Merged.0v24.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 1DPF - Log2FC',col='red',main='Watersnake',xlim=c(-6,6),ylim=c(-6,6))
# points(Ner.Merged.0v24.FC.notsig$RNA_Log2FC,Ner.Merged.0v24.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.2673,a=-0.1977)
# text(x=-5,y=4.5,labels = expression(R^2 == 0.08717))
# text(x=-5,y=4,labels='p < 0.05')
# abline(b=0.27566,a=-0.55208 ,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.02937),col='blue')
# text(x=-5,y=3,labels='p < 0.001',col='blue')
# 
# 
# 
# #python b
# plot(Py.Merged.0v24.0v96.FC.sig$RNA_Log2FC,Py.Merged.0v24.0v96.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 4DPF - Log2FC',col='red',xlim=c(-6,6),ylim=c(-6,6))
# points(Py.Merged.0v24.0v96.FC.notsig$RNA_Log2FC,Py.Merged.0v24.0v96.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.40157,a=-0.78841)
# text(x=-5,y=4.5,labels = expression(R^2 == 0.09673))
# text(x=-5,y=4,labels='p < 0.001')
# abline(b=0.48128,a=-1.09599,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.07928),col='blue')
# text(x=-5,y=3,labels='p < 0.001',col='blue')
# 
# #rattlesnake b
# plot(CV.Merged.0v24.0v96.FC.sig$RNA_Log2FC,CV.Merged.0v24.0v96.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 4DPF - Log2FC',col='red',xlim=c(-6,6),ylim=c(-6,6))
# points(CV.Merged.0v24.0v96.FC.notsig$RNA_Log2FC,CV.Merged.0v24.0v96.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.1883,a=-0.3589)
# text(x=-5,y=4.5,labels = expression(R^2 == 0.01122))
# text(x=-5,y=4,labels='p = 0.1279')
# abline(b=0.16382,a=-0.44073,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.005263),col='blue')
# text(x=-5,y=3,labels='p = 0.08861',col='blue')
# 
# 
# #watersnake b
# plot(Ner.Merged.0v24.0v96.FC.sig$RNA_Log2FC,Ner.Merged.0v24.0v96.FC.sig$Prot_Log2FC,xlab='RNA Fasted vs 1DPF - Log2FC',ylab='Protein Fasted vs 4DPF - Log2FC',col='red',xlim=c(-6,6),ylim=c(-6,6))
# points(Ner.Merged.0v24.0v96.FC.notsig$RNA_Log2FC,Ner.Merged.0v24.0v96.FC.notsig$Prot_Log2FC,col='black')
# abline(b=0.06816,a=-0.54930)
# text(x=-5,y=4.5,labels = expression(R^2 == 0))
# text(x=-5,y=4,labels='p = 0.6019')
# abline(b=0.12910,a=-1.11254,col='blue')
# text(x=-5,y=3.5,labels = expression(R^2 == 0.00227),col='blue')
# text(x=-5,y=3,labels='p = 0.1441',col='blue')
# par(mfrow=c(1,1))
# 
# 
# 
# # ###
# # ### Experimenting w/ RNA-prot correlation within individuals
# # ###
# # 
# # Py.1.normPSM <- Py.normPSM[,1:2]
# # Py.1.normRNA <- Py.RNA.NormCounts[,1:2]
# # Py.1.merged <- merge(Py.1.normPSM,Py.1.normRNA,all.x=T)
# # Py.1.merged[,2:3] <- log2(Py.1.merged[,2:3])
# # plot(Py.1.merged$AI6_0h,Py.1.merged$py1_X..PSMs)
# # Py.1.merged <- Py.1.merged[which(!is.infinite(Py.1.merged$py1_X..PSMs)),]
# # Py.1.merged <- Py.1.merged[which(!is.infinite(Py.1.merged$AI6_0h)),]
# # 
# # Py.1.lm <- lm(formula = py1_X..PSMs ~ AI6_0h, data = Py.1.merged)
# # summary(Py.1.lm)
# # abline(b=0.1337,a=1.4539)
# # text(x=1,y=6,labels = expression(R^2 == 0.05461))
# # text(x=1,y=5.5,labels='p < 0.001')
# # 
# # 
# # Py.2.normPSM <- Py.normPSM[,c(1,3)]
# # Py.2.normRNA <- Py.RNA.NormCounts[,c(1,3)]
# # Py.2.merged <- merge(Py.2.normPSM,Py.2.normRNA,all.x=T)
# # Py.2.merged[,2:3] <- log2(Py.2.merged[,2:3])
# # plot(Py.2.merged[,3],Py.2.merged[,2])
# # 
# # Py.2.merged <- Py.2.merged[which(!is.infinite(Py.2.merged$py2_X..PSMs)),]
# # Py.2.merged <- Py.2.merged[which(!is.infinite(Py.2.merged$AJ6_0h)),]
# # 
# # Py.2.lm <- lm(formula = py2_X..PSMs ~ AJ6_0h, data = Py.2.merged)
# # summary(Py.2.lm)
# # abline(b=0.07499,a=1.81696)
# # text(x=1,y=6,labels = expression(R^2 == 0.0183))
# # text(x=1,y=5.5,labels='p < 0.001')
# # 
# # Py.3.normPSM <- Py.normPSM[,c(1,4)]
# # Py.3.normRNA <- Py.RNA.NormCounts[,c(1,4)]
# # Py.3.merged <- merge(Py.3.normPSM,Py.3.normRNA,all.x=T)
# # Py.3.merged[,2:3] <- log2(Py.3.merged[,2:3])
# # plot(Py.3.merged[,3],Py.3.merged[,2])
# # 
# # Py.4.normPSM <- Py.normPSM[,c(1,5)]
# # Py.4.normRNA <- Py.RNA.NormCounts[,c(1,5)]
# # Py.4.merged <- merge(Py.4.normPSM,Py.4.normRNA,all.x=T)
# # Py.4.merged[,2:3] <- log2(Py.4.merged[,2:3])
# # plot(Py.4.merged[,3],Py.4.merged[,2])
# # 
# # 
# # Py.5.normPSM <- Py.normPSM[,c(1,6)]
# # Py.5.normRNA <- Py.RNA.NormCounts[,c(1,6)]
# # Py.5.merged <- merge(Py.5.normPSM,Py.5.normRNA,all.x=T)
# # Py.5.merged[,2:3] <- log2(Py.5.merged[,2:3])
# # plot(Py.5.merged[,3],Py.5.merged[,2])
# # Py.5.merged <- Py.5.merged[which(!is.infinite(Py.5.merged$py5_X..PSMs)),]
# # Py.5.merged <- Py.5.merged[which(!is.infinite(Py.5.merged$AF2_0h)),]
# # 
# # Py.5.lm <- lm(formula = py5_X..PSMs ~ AF2_0h, data = Py.5.merged)
# # summary(Py.5.lm)
# # abline(b=0.1219,a=1.3325)
# # text(x=1,y=6,labels = expression(R^2 == 0.04755))
# # text(x=1,y=5.5,labels='p < 0.001')
# # 
# # 
# # Py.7.normPSM <- Py.normPSM[,c(1,7)]
# # Py.7.normRNA <- Py.RNA.NormCounts[,c(1,7)]
# # Py.7.merged <- merge(Py.7.normPSM,Py.7.normRNA,all.x=T)
# # Py.7.merged[,2:3] <- log2(Py.7.merged[,2:3])
# # plot(Py.7.merged[,3],Py.7.merged[,2])
# # 
# # Py.7.normPSM <- Py.normPSM[,c(1,8)]
# # Py.7.normRNA <- Py.RNA.NormCounts[,c(1,8)]
# # Py.7.merged <- merge(Py.7.normPSM,Py.7.normRNA,all.x=T)
# # Py.7.merged[,2:3] <- log2(Py.7.merged[,2:3])
# # plot(Py.7.merged[,3],Py.7.merged[,2])
# # 
# # Py.8.normPSM <- Py.normPSM[,c(1,9)]
# # Py.8.normRNA <- Py.RNA.NormCounts[,c(1,9)]
# # Py.8.merged <- merge(Py.8.normPSM,Py.8.normRNA,all.x=T)
# # Py.8.merged[,2:3] <- log2(Py.8.merged[,2:3])
# # plot(Py.8.merged[,3],Py.8.merged[,2])
