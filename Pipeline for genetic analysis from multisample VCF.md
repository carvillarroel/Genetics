# Pipeline for genetic analysis from multisample VCF #
## Working with VCF ##

First vcftools is used to filter SNPs, individuals, etc.
For the first analyisis we will use SNPs that were called in all individuals, and we only retain biallelic sites  

```bash
conda install vcftools perl-vcftools-vcf
vcftools --max-missing 1 --max-alleles 2 --vcf FILE.vcf --recode --recode-INFO-all --out OUTFILE
```


## Importing VCF into R ##

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
This phylogenetic tree can be viewed and modified using iTol website for example. However we will generate a more robust phylogenetic tree using Maximum Likehood method implemented in the program IQTREE.

```bash
conda install iqtree
git clone https://github.com/edgardomortiz/vcf2phylip.git
#First we transform the VCF file to phylip format using a custom python script developed by edgardomortiz, and then we use it for IQTREE
python vcf2phylip -i OUTFILE.recode.vcf
iqtree -s OUTFILE.min4.recode.min4.phy -st DNA -o 1105.1_Nahuelbuta -m GTR+ASC -nt 8 #outgroup is 1105.1_Nahuelbuta as is a different species
iqtree -s OUTFILE.min4.recode.min4.phy.varsites.phy -st DNA -o 1105.1_Nahuelbuta -m GTR+ASC -nt 8 -bb 1000 -redo
```
The file OUTFILE.min4.recode.min4.phy.varsites.phy.contree is the ML phylogenetic tree that can also be viewed in iTol.









