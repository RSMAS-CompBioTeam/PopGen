# PopGen
Scripts for CompBio-athon PopGen instructional sessions: SNP analyses, Population Differentiation analyses, Outlier analyses &amp; GWAS, eQTL analyses

In this tutorial we will learn:
  1) How to call SNPs
  2) Filtering using vcftools
  3) Plotting a PCA
  4) Identifying FST outliers


# How to call SNPs
Once you have aligned your reads to a reference genome or transcriptome, there are many different programs that can be used to identify Single Nucleotide Polymorphisms (SNPs). Some of the most widely used include:
  * GATK (https://software.broadinstitute.org/gatk/)
  * Freebayes (https://github.com/ekg/freebayes)
  * Samtools mpileup (http://samtools.sourceforge.net/mpileup.shtml)
  * ANGSD (for genotype likelihoods) (https://github.com/ANGSD/angsd)

Today, we will be using one of the more simple pipelines - **samtools**. We start with a set of bam alignment files and a fasta file of the reference genome. We'll use a set downloaded from dropbox along with some scripts we will be using:

```bash
wget --no-check-certificate https://www.dropbox.com/sh/9l812bwwnorwfom/AAAF3VYK1VmVLMlXzeR-aWtFa?dl=1
unzip file
```

The script `snp.sh` shows how we can use samtools to identify SNPs in our dataset. Let's execute it:

```bash
bsub < snp.sh
```

This usually takes a while - up to a few days depending on how much data you have. The output will be a file in Variant Call Format (VCF). Go ahead and look at this file (here is a reference for what the fields mean: https://samtools.github.io/hts-specs/VCFv4.2.pdf). This VCF will serve as the raw data for all of our downstread SNP analyses.

# Filtering using vcftools
The raw VCF file contains a lot of information. We have SNPs, but we also have INDELS. We have variants that are only seen in one individual and we have regions of the genome where we don't even have coverage for many individuals. We have variants of very low quality. The next step is to filter out the variants that we don't want. the best tool for this is **vcftools** (http://vcftools.sourceforge.net/). We'll do this in interactive mode on Pegasus:
```bash
bsub -Is -q interactive bash
```
Pegasus already has vcftools, so we can just load the module. We also have to load the correct perl module so that it works.
```bash
module unload perl
module load perl/5.22.1
module load vcftools/0.1.13
```
Let's just execute the most simple commange and see what we get:
```bash
vcftools --vcf raw_snps.vcf
```
This tells you how many individuals and how many variants there are in your vcf. It does not save anything and it does not do any filtering. Let's try to do a simple filtering step:
```bash
vcftools --vcf raw_snps.vcf --max-missing 0.5
```
This means that we have to have genotypes in at least 50% of our individuals. How many SNPs do you have now? Notice that we still didn't output a new file. To do that:
```bash
vcftools --vcf raw_snps.vcf --recode --out WIFL_filtered
```
There are many different filters you can work with. Consult the VCFtools manual for the full list. Here are a few of my favorites:

  --remove-indels | removes INDELS so only SNPs remain
  
  --min-alleles | minimum numbers of allelic states observed
  
  --max-alleles | maximum number of allelis states
  
  --maf | minor allele frequency cutoff
  
  --minGQ | minimum genotype quality (phred score)
  
  --minDP | minimum read depth
  
In order to save the vcf, you must use the following options:

  --recode
  
  --out filtered

Filter this set of SNPs however you want, we'll use them to make a PCA and for outlier analysis. You can even save multiple filters is you want to compare the downstream analysis!

# Analyzing SNP data
Now that we've filtered our SNP data, we are ready to carry out downstream statistical analyses, which we will do in RStudio. There are a couple of R scripts in this Github repository that we'll use for these analyses. An easy way for you to get these scripts onto your computer (which allows us to show you more cool ways to use Github!) is to fork and clone this repository. Creating a fork of this repository makes a copy of it that is owned by you. Once you make this forked copy, you can then clone it (i.e. download it) to your local computer and work with the files inside. Because this copy of the repository is owned by you, you can also make changes and then push them back to Github. 

First, click the 'Fork' button at the top right of this page, and select your Github account. Next, open RStudio and click File > New Project. This time, select 'Version Control'. Paste in the web address to your new forked repository, call it 'PopGen', and create it within your `~/github/` directory.

Cool - now you have a new Rproject that contains the scripts we are going to use in the next part of the analysis! Git is tracking everything already so you can keep track of your analysis as it develops from here.


# Principal Components Analysis
One of the most basic ways to visualize mulitvariate data is with a principal components analysis (PCA). Here we'll make one based on the filtered SNP file you made in the last step. Since we'll do this in R on our local computers, we first need to get our vcf files from Pegasus into our local repository (make sure to change username and paths):
```bash
scp username@pegasus.ccs.miami.edu:/path/to/vcf ~/github/PopGen/
```
Now just step through the `SNPRelate.R` script to make your PCA!
Note: SNPRelate is a great package that can do many other things besides PCA visualization, including calculating FST and relatedness. If you have extra time, you can play around with these functions. Here is a great tutorial: https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html

If you still have extra time, you can play around with a spatial structure analysis using TESS3 (https://github.com/cayek/TESS3/tree/master/tess3r). Install the R package and follow the tutorial for help with the code.

# Outlier Analysis
There are many software programs designed to do outlier analysis. We'll use OutFLANK because it is easy to implement and its non-parametric approach is good for organisms with complex demography. Work through the `OutFLANK.R` script step by step.
