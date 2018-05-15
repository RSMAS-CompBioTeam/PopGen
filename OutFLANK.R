library(OutFLANK)
library(vcfR)

###Paths to data
vcf.path = "~/Desktop/compbio/filtered.recode.vcf"
meta.path = "~/Desktop/compbio/WIFL_meta.txt"

###VCF to OutFLANK format - info on this found in OutFLANK vignette
data <- read.vcfR(vcf.path)
geno <- extract.gt(data)
position <- getPOS(data)
chromosome <- getCHROM(data)
position <- 1:nrow(geno)
G <- matrix(NA,nrow=nrow(geno),ncol=ncol(geno))
G[geno %in% c("0/0", "0|0")] <- 0
G[geno  %in% c("0/1", "1/0", "1|0", "0|1")] <- 1
G[geno %in% c("1/1", "1|1")] <- 2
G[is.na(G)] <- 9
table(as.vector(G))

###population data
meta <- read.delim(meta.path)
pops <- meta$Group[match(colnames(geno),meta$Sample)]

###Make FST matrix
fst <- MakeDiploidFSTMat(t(G),locusNames=position,popNames=pops) # This takes a couple minutes
plot(fst$FSTNoCorr,fst$FST) # This is a check to see if we can fit Chi-Sq -should look like a straight line

###Run OutFLANK
OF <- OutFLANK(fst,LeftTrimFraction=0.05,RightTrimFraction=0.05,
         Hmin=0.1,NumberOfSamples=164,qthreshold=0.1)
OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=F,RightZoomFraction=0.01,titletext=NULL) #Plot distribution
OutFLANKResultsPlotter(OF,withOutliers=T,
                       NoCorr=T,Hmin=0.1,binwidth=0.005,
                       Zoom=T,RightZoomFraction=0.01,titletext=NULL) #Zoom in on right tail

###Outliers
P1 <- pOutlierFinderChiSqNoCorr(fst,Fstbar=OF$FSTNoCorrbar,
                                dfInferred=OF$dfInferred,qthreshold=0.1,Hmin=0.1)
outliers <- P1$OutlierFlag==TRUE
table(outliers)
plot(P1$He,P1$FST,col=rgb(0,0,0,alpha=0.1))
points(P1$He[outliers],P1$FST[outliers],col="magenta")

###Manhattan Plot
plot(P1$LocusName,P1$FST,xlab="Position",ylab="FST",col=rgb(0,0,0,alpha=0.1))
points(P1$LocusName[outliers],P1$FST[outliers],col="magenta")
