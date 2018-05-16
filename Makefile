all: myPCA.png FST.png

myPCA.png: SNPRelate.R WIFL_filtered.vcf
	R --vanilla < SNPRelate.R
	
FST.png: OutFLANK.R WIFL_filtered.vcf
	R --vanilla < OutFLANK.R

# WIFL_filtered.vcf: WIFL.vcf
#	vcftools --vcf raw_snps.vcf --max-missing 0.5 --recode --out WIFL_filtered