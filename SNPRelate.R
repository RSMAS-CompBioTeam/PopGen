library(SNPRelate)

### Relative paths to data within the PopGen repository
vcf.path = "filtered.recode.vcf"
meta.path = "WIFL_meta.txt"

###Convert to gds file
snpgdsVCF2GDS(vcf.path,"genos.gds",method="biallelic.only") # This take a while
snpgdsSummary("genos.gds")

###Open gds file and read in population data
genofile <- snpgdsOpen("genos.gds")
meta <- read.delim(meta.path)


###LD-based SNP pruning
set.seed(1000)
snpset <- snpgdsLDpruning(genofile,ld.threshold=0.5,autosome.only=F) # You can change LD threshold
snpset.id <- unlist(snpset)

###Calculate PCA
pca <- snpgdsPCA(genofile,snp.id=snpset.id,autosome.only=F)
pc.percent <- pca$varprop*100
head(round(pc.percent,2)) # Proportion of variance explained by each PC axis

###Make dataframe and plot
tab <- data.frame(sample.id=pca$sample.id,
                  pop=factor(meta$Group)[match(pca$sample.id,meta$Sample)],
                  EV1=pca$eigenvect[,1],
                  EV2=pca$eigenvect[,2],
                  stringsAsFactors=F)
plot(tab$EV2,tab$EV1,col=as.integer(tab$pop),xlab="eigenvector 2",ylab="eigenvector 1")
legend("topright",legend=levels(tab$pop),pch="o",col=1:nlevels(tab$pop))

snpgdsClose(genofile)

### Save plot as a png image
#png(filename = "myPCA.png", width = 4, height = 4, units = "in", res = 300)
#plot(tab$EV2,tab$EV1,col=as.integer(tab$pop),xlab="eigenvector 2",ylab="eigenvector 1")
#legend("topright",legend=levels(tab$pop),pch="o",col=1:nlevels(tab$pop))
#dev.off()
