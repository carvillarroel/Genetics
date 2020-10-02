# Genetic analyses from whole genome resequencing #
## Phylogenetic tree and population structure
### Working with VCF ###

First vcftools is used to filter SNPs, individuals, etc.
For the first analyisis we will use SNPs that were called in all individuals, and we only retain biallelic sites  

```bash
conda install vcftools perl-vcftools-vcf
vcftools --max-missing 1 --max-alleles 2 --vcf FILE.vcf --recode --recode-INFO-all --out OUTFILE
```


### Importing VCF into R ###

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
This phylogenetic tree can be viewed and modified using iTol website for example. However we will generate a more robust phylogenetic tree using Maximum Likehood method implemented in the linux program IQTREE. First we transform the VCF file to phylip format using a custom python script developed by edgardomortiz, [vcf2phylip](https://github.com/edgardomortiz/vcf2phylip), and then we use it for IQTREE

```bash
conda install -c bioconda iqtree
git clone https://github.com/edgardomortiz/vcf2phylip.git
#First we transform the VCF file to phylip format using a custom python script developed by edgardomortiz, and then we use it for IQTREE
python vcf2phylip -i OUTFILE.recode.vcf
iqtree -s OUTFILE.min4.recode.min4.phy -st DNA -o 1105.1_Nahuelbuta -m GTR+ASC -nt 8 #outgroup is 1105.1_Nahuelbuta as is a different species
iqtree -s OUTFILE.min4.recode.min4.phy.varsites.phy -st DNA -o 1105.1_Nahuelbuta -m GTR+ASC -nt 8 -bb 1000 -redo
```
The file OUTFILE.min4.recode.min4.phy.varsites.phy.contree is the ML phylogenetic tree that can also be viewed in iTol.

### Exploring population structure ###
We use STRUCTURE to estimate the number of populations ocurring in our dataset. Structure is recommended to work with a subset of SNPs (~10k snps), which we can obtain by filtering out SNPs that are in linkage disequilibrium (LD) using the software PLINK.
As for determining population structure we do not need an outgroup, we can remove it using vcftools (preferentially using the first unfiltered vcf).

```bash
Create a text file with the name of the individual to remove (remove.txt)
vcftools --remove remove.txt --vcf FILE.vcf --recode --recode-INFO-all --non-ref-ac-any 1 --out onlyeub
```

Now to use PLINK, we modify the VCF file to write names to each SNP (requirement if we want to use the list that PLINK produces with the SNPs that are in LD)
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

Now we use this file to retain these positions from the vcf with vcftools. However we need to change the underscore character that is located prior to the SNP position to a tab character
For example: CBS12357_Chr12_polished_838562 --> CBS12357_Chr12_polished 838562

```bash
sed 's/polished_/polished\t/' onlyeub_ldfilter.prune.in > onlyeub_ldfilter.prune.in.vcftools
vcftools --vcf onlyeub.recode.vcf --positions onlyeub_ldfilter.prune.in.vcftools --recode-INFO-all --recode --out onlyeub_ldfilter
```
If we still have more than 20k SNPs, we can filter SNPs by distance using vcftools --thin option. Here we will filter out SNPs that are closer than 500 bp

```bash
vcftools --vcf onlyeub_ldfilter.recode.vcf --thin 500 --recode --recode-INFO-all --out onlyeub_ldfilter_thinned
```

Now to use this VCF file in structure we need to change its format. Here we use STACKS to convert the VCF file to structure format.
```bash
conda install structure stacks
mkdir onlyeub_structure
populations -V onlyeub_ldfilter_thinned.recode.vcf -O onlyeub_structure --structure
cd onlyeub_structure
tail -n +2 onlyeub_ldfilter_thinned.recode.p.structure > onlyeub_structure.txt
```
The tail command is just to delete the first line of the file we need. Also stacks (populations) will generate a few files that we don't need.
Besides the input, Structure needs two files with the options for the run: a mainparams and extraparams files. Extraparams can be empty, but we need to populate the mainparams file

```bash
touch extraparams
touch mainparams
nano mainparams
```

The mainparams should look like this:
```bash
#define OUTFILE k2
#define INFILE onlyeub_structure.txt
#define NUMINDS 107
#define NUMLOCI 15003
#define LABEL 1 
#define POPDATA 0 
#define POPFLAG 0 
#define LOCDATA 0 
#define PHENOTYPE 0 
#define MARKERNAMES 1 
#define MAPDISTANCES 0 
#define ONEROWPERIND 0 
#define PHASEINFO 0 
#define PHASED 0 
#define RECESSIVEALLELES 0 
#define EXTRACOLS 1
#define MISSING 
#define PLOIDY 2
#define MAXPOPS 2
#define BURNIN 10000
#define NUMREPS 100000


#define NOADMIX 0
#define LINKAGE 0
#define USEPOPINFO 0

#define LOCPRIOR 0
#define INFERALPHA 1
#define ALPHA 1.0
#define POPALPHAS 0 
#define UNIFPRIORALPHA 1 
#define ALPHAMAX 10.0
#define ALPHAPROPSD 0.025


#define FREQSCORR 1 
#define ONEFST 0
#define FPRIORMEAN 0.01
#define FPRIORSD 0.05


#define INFERLAMBDA 0 
#define LAMBDA 1.0
#define COMPUTEPROB 1 
#define PFROMPOPFLAGONLY 0 
#define ANCESTDIST 0 
#define STARTATPOPINFO 0 
#define METROFREQ 10


#define UPDATEFREQ 1 
```

We mainly need to change OUTFILE,INFILE,NUMINDS,NUMLOCI,MAXPOPS,BURNIN,NUMREPS to reflect our file and desired setup. Importantly, MAXPOPS is the number of populations to test, and we want to explore many values of MAXPOPS to later estimate the best number of populations that can explain our genetic data. Fortunately we do not need make a mainparams file for each different iteration of Structure we want to run, as the command line version can accept options that will override the mainparams file (but still needs this file to run!).
Options that can be added in command line:
```
-m mainparams file
-e extraparams file 
-K MAXPOPS = number of populations to test
-L NUMLOCI = total number of SNPs in the vcf file used for structure input
-N NUMINDS = total number of individuals in the file
-i input file 
-o output file
```
So first we check that everything is ok by doing a test run  (change the parameters to match your file!):

```bash
structure -m mainparams -p extraparams -k 2 -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_res_K2_rep1.txt
```

If structure runs without problems, we can cancel the test run to run structure in a loop to generate multiple replicate runs for each value of K that we want to perform. Importantly we need to give some waiting time between runs as structure produces a SEED value from the system clock, and we dont want to use the same seed for the different replicates.

```bash
for m in {2..7}
   do
   structure -m mainparams -p extraparams -K $m -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_K.${m}.rep1 &> log_K.${m}.rep1.txt&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_K.${m}.rep2 &> log_K.${m}.rep2.txt&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_K.${m}.rep3 &> log_K.${m}.rep3.txt&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_K.${m}.rep4 &> log_K.${m}.rep4.txt&
   sleep 10
   structure -m mainparams -p extraparams -K $m -L 15003 -N 107 -i onlyeub_structure.txt -o onlyeub_structure_K.${m}.rep5 &> log_K.${m}.rep5.txt
   sleep 10  
   done 
done
```

The above loop will run structure from K=2 to K=7, doing 5 replicates for each K in parallel. As the value of K increases, also does the time to run structure.   


