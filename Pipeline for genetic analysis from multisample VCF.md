### Pipeline for genetic analysis from multisample VCF ###
# Working with VCF #

First vcftools is used to filter SNPs, individuals, etc.
For the first analyisis we will use SNPs that were called in all individuals, and only biallelic sites  

```bash
conda install vcftools perl-vcftools-vcf
vcftools --max-missing 1 --max-alleles 2 --vcf FILE.vcf --recode --recode-INFO-all --out OUTFILE
```


# Importing VCF into R #

In R we will use different packages to perform an array of genetic analysis. First we will import the VCF with the package vcfR, and then transform it to a genlight object (a simplified version of the VCF with only allele information as binary code)

```r
library(vcfR)
library(adegenet)
vcf<-read.vcfR("OUTFILE.recode.vcf")
x <- vcfR2genlight(vcf)
ploidy(x) <- 2 
x@ind.names
write.csv(as.matrix(x@ind.names),file="indivuals.csv") ## Export individuals names to complete with metainformation (populations for example)
```

Now this genlight object can be processed with adegenet, poppr, and other genetic analysis packages. First we will perform a phylogenetic tree using a bitwise distance matrix built from SNP data
```r
library(poppr)
tree<-aboot(x,tree = "nj", distance = bitwise.dist, sample = 1000, showtree = F, cutoff = 50, quiet = T)

write.nexus(tree,file="tree.nexus")
```
This phylogenetic tree can be viewed and modified in the using iTol website. However we will generate a more robust phylogenetic tree later.


