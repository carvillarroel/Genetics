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

## Exploring population structure ##
First we use STRUCTURE to obtain populations from our dataset. Structure is recommended to work with a subset of SNPs (aprox 10k), which we can first obtain by filtering SNPs that are in linkage disequilibrium (LD) using the software PLINK.
At this step we do not need an outgroup, so we can remove it using vcftools (preferentially using the unfiltered vcf)
Create a text file with the name of the individual to remove (remove.txt)
vcftools --remove remove.txt --vcf FILE.vcf --recode --recode-INFO-all --non-ref-ac-any 1 --out onlyeub

Now to use PLINK, we first need to modify the VCF to give names to each SNP (requirement if we want to use the list that PLINK produces with the SNPs that are in LD)
We will use the following script available in gist:
https://gist.github.com/janxkoci/25d495e6cb9f21d5ee4af3005fb3c77a#file-plink_pruning_prep-sh

```bash
wget https://gist.github.com/janxkoci/25d495e6cb9f21d5ee4af3005fb3c77a/raw/a18b0d37f6bbf0a354bd928b158efb0bd85bd916/plink_pruning_prep.sh
conda install tabix bcftools plink
./plink_pruning_prep.sh onlyeub.recode.vcf
plink --vcf onlyeub.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.2 --out onlyeub_ldfilter

#Check if the prune.in file contain SNP ids#
head onlyeub_ldfilter.prune.in
```

Now we use this file to retain these positions from the vcf with vcftools. However we need to change the "_" character prior to the SNP position to a TAB, to obtain this format: CHROMOSOME<TAB>POSITION
CBS12357_Chr12_polished_838562 to CBS12357_Chr12_polished 838562

```bash
sed 's/polished_/polished\t/' onlyeub_ldfilter.prune.in > onlyeub_ldfilter.prune.in.vcftools
vcftools --vcf onlyeub.recode.vcf --positions onlyeub_ldfilter.prune.in.vcftools --recode-INFO-all --out onlyeub_ldfilter
```
If we still have more than 20k SNPs, we can filter SNPs by distance using vcftools --thin option. Here we will filter out SNPs that are closer than 500 bp

```bash
vcftools --vcf onlyeub_ldfilter.recode.vcf --thin 500 --recode --recode-INFO-all --out onlyeub_ldfilter_thinned
```

Now to use this VCF file in structure we need to change its format. Here we use STACKS to obtain to convert the VCF file.
```bash
conda install structure stacks
mkdir onlyeub_structure
populations -V onlyeub_ldfilter_thinned.recode.vcf -O onlyeub_structure --structure
cd onlyeub_structure
tail -n +2 onlyeub_ldfilter_thinned.recode.p.structure > onlyeub_structure.txt
```
The tail command is just to delete the first line of the file we need. Also stacks (populations) will generate a few files that we don't need.
Besides the input, Structure needs two files with the options for the run: a mainparams and extraparams files.






