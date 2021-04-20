# Genetic analyses from whole genome resequencing #
## Phylogenomics analysis with TREEMIX

The program TREEMIX allows us to examine possible events of gene flow that happened a long time ago between the ancestral populations of TODAY sequenced individuals.
Install TREEMIX using conda

To perform TREEMIX analysis we will start from the first unfiltered vcf we created and modify it with a series of steps (many of which already used when we created the dataset for STRUCTURE)

As we want to look at old gene flow events, we are not interested in recent gene flow events, then we exclude from this analysis admixed individuals, and we keep preferably "pure" individuals
1) We will create a tabulated file with the individuals (column 1) that we keep and their corresponding subpopulation (column 2)(ej: PB1-Karukinka), and name it popstacks.txt

2) Exclude individuals not analyzed from the vcf using vcftools:

```bash
vcftools --vcf filtered_snps.vcf --keep keep.txt --recode-INFO-all --recode --max-alleles 2 --non-ref-ac-any 1 --max-missing 1 --out treemix_vcf_f1
```

3) Filter sites in LD using the script prunning_ld.sh and PLINK (As we did for STRUCTURE)

```bash

./plink_pruning_prep.sh treemix_vcf_f1.recode.vcf
plink --vcf treemix_vcf_f1.recode_annot.vcf --double-id --allow-extra-chr --indep-pairwise 50 10 0.2 --out treemix_ldfilter
sed 's/polished_/polished\t/' treemix_ldfilter.prune.in > treemix_ldfilter.prune.in.vcftools
vcftools --vcf treemix_vcf_f1.recode.vcf --positions treemix_ldfilter.prune.in.vcftools --recode-INFO-all --recode --out treemix_vcf_f2

```
4) Use STACKS to create a TREEMIX file from the filtered VCF. Here we use the popstacks.txt file we created before

```bash
populations -V treemix_vcf_f2.recode.vcf -O treemix -M popstacks2.txt --treemix
```
5) Delete the first line of the STACKS output, and gzip it
```bash
sed -i '1d' treemix_vcf_f2.recode.p.treemix
gzip treemix_vcf_f2.recode.p.treemix
```
6) We need to give TREEMIX an hypothesis of how many migration (gene flow) events ocurred in our datasets, but if we dont know a proper number we could test for the "best" number of migration events, similar to the best K of STRUCTURE.
To do so we will run TREEMIX in a loop to test migrations from for example m=0 to m=10. Create the following bash file:

```bash
nano treemix_batch.sh


---- in NANO --------
for m in {0..10}
   do
   for i in {1..10}
      do
      treemix \
         -i $1 \
         -o all_bestm.${i}.${m} \
         -m ${m} \
         -k 500 \
         -noss \
        
      done 
done

------------------
```

./treemix_batch.sh treemix.

7) Now we use the R package OPTM to find the best values of m (https://cran.r-project.org/web/packages/OptM/index.html)

```R

library("optM")
test.optM = optM("treemix/")  #HERE THE FOLDER WITH ALL OUTPUTS FROM TREEMIX 
plot_optM(test.optM, method = "Evanno")
```
8) With the best m we now perform several runs of TREEMIX which will allow us to add boostrap information using the R package BITE (https://github.com/marcomilanesi/BITE)

Download BITE from github and decompress. Then decompress the treemix_scripts.tar.bz2 file inside the inst/script folder. We will use the script Treemix_bootstrap.sh
For this script to work, we need phylip installed (conda install phylip), and we need to give the scrip the location of the binary file "consense" after we installed with CONDA.

Then we run the script with the arguments in order:
bash Treemix_boostrap.sh TREEMIX_FILE MIGRATIONS CPUs BINs OUTGROUP_POPULATION BOOSTRAPS CONSENSE_LOCATION OUTPUT_NAME


```bash
bash Treemix_bootstrap.sh ../../../../treemix_vcf_f2.recode.p.treemix.gz 1 4 500 UVA 100 ~/miniconda3/bin/consense ite1_1mig500bin
```

Then in R, install BITE from source (using the decompressed files), but before its dependencies (check website)

```R

library(BITE)

pdf("Test_Treemix_bootstrap.pdf", 10, 7)
treemix.bootstrap(in.file = "ite1_1mig500bin", 
	  out.file = "test", 
    phylip.file = "ite1_1mig500bin_outtree.newick", 
    nboot = 100, 
    cex=0.5, xmin = -0.005, disp = 0.001, 
    boot.legend.location = "top", xbar = 0.05)
dev.off()

}
```
