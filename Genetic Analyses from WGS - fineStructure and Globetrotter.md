# Chromopainter - fineSTRUCTURE - GLOBETHROTTER pipeline

To run these three programs we need a phased SNP file.  
For this the phasing pipeline included in the program GERMLINE (1.5.3) was used, which uses BEAGLE version 3.0.4 (there are several new versions of BEAGLE but the output of these versions wont work with chromopainter) and PLINK to perform haplotype phasing

Programs:
```bash
conda install snpsift plink
wget https://people.maths.bris.ac.uk/~madjl/finestructure/fs_4.1.1.zip
wget https://github.com/gusevlab/germline/blob/master/phasing_pipeline.tar.gz
wget https://faculty.washington.edu/browning/beagle/recent.versions/beagle_3.0.4_05May09.zip <- JUST NEED TO COPY BEAGLE.JAR TO THE PHASING PIPELINE FOLDER
```
Steps

1.- Split VCF file (105 eubayanus strains, no uvarum) by chromosome using Sift
```bash
SnpSift split alleub_f1.recode.vcf
```
2.- Convert VCF files to PLINK .ped / .map format
```bash
plink --vcf alleub_f1.recode.CBS12357_Chr01_polished.vcf --recode12 --allow-extra-chr --double-id --geno 1 --out chr1
```
iterate for each chromosome

PHASING (NOT NECESSARY FOR HAPLOIDS!)
3.- Run script run.sh from the phasing_pipeline package for each PLINK file (IT NEEDS BEAGLE.jar IN THE RUN.SH FOLDER)
```bash

for i in {1..16};do bash run.sh chr${i}.ped chr{i}.map chr{i};done
```
4.- Convert phased files to chromopainter format
```bash

for i in {1..16};do perl plink2chromopainter.pl -p=chr${i}.phased.ped -m=chr{i}.phased.map -o=chr{i}.chromopainter -f;done

```

##FOR HAPLOIDS
USE THIS SCRIPT TO CHANGE CHROMOPAINTER INPUT TO HAPLOID (Example here Lachancea with chromosomes names)
```bash
for i in  LACI0A LACI0B LACI0C LACI0D LACI0E LACI0F LACI0G LACI0H;do (awk 'NR == 1  { print $1 /2 }' chr_${i}.chromopainter; sed '2,3!d' chr_${i}.chromopainter;sed '1,3d' chr_${i}.chromopainter|sed  '0~2d')| cat > chr_${i}.chromopainter.haploid;done
```


5.- Create recombination file for each chromosome (edit line 47 in the perl script makeuniformrecfile.pl to give a constant value of 4/1,000,000 (0.000004)) which means that in all the following calculations we will use a constant recombination rate of 0.4 cM/Kb (which is the average in S. cerevisiae (Cubillos et al, 2011)). Iterate this perl script over the files created in step 5
```bash

	for i in {1..16};do perl makeuniformrecfile.pl  chr{i}  chr{i}_rec;done
```
6.-	Run chromopainter (v2) in mode ALLvsALL (all individuals will be painted with each other individual)
	To run 4 processes:
		for i in {1..2};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {3..4};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {5..6};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {7..8};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {9..10};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {11..12};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {13..14};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &
		for i in {15..16};do ../../fs_4.0.1/fs_linux_glibc2.3 chromopainter -g ../chr${i}.chromopainter -r ../chr${i}_rec -t idfile.txt -o cp_chr${i} -a 0 0;done &

		Obs1: idfile is a tabulated file with the name of each strain and their corresponding population, and whether to include this individual in the calculation. Should have the same order that the VCF and PLINK files. e.g:
		CL1001.1<TAB>PB1<TAB>1

		Obs2: In this step it will be better to change the name of the individuals if they start with a number (eg. 1001.1_TDP change it to CL1001.1), as finestructure might give problems with names that start with numbers later
		Obs3: option -a 0 0 means ALL vs ALL mode
		Obs4: FOR HAPLOIDS USE ADDITIONAL OPTION -j 1


7.-	Run chromocombine to combine all chromopainter per chromosome outputs into a combined genome output (put all of them in a folder e.g. datos)
		../../fs_4.0.1/fs_linux_glibc2.3 chromocombine -d datos

8.-	Run finestructure on the combined dataset 
		../../fs_4.0.1/fs_linux_glibc2.3 finestructure  -x 100000 -y 100000 -z 1000 output.chunkcounts.out out.105eubs.mcmc.xml
		Visualize result with finestructure GUI or...
		Use finestructure R scripts:
		Obs: To plot the results with finestructure R scripts, I had to rename the individuals to this format:
		eg:	CL1001.1 to PBI1001 or CL210.1 to PBIII210
		This is because the script developed by finestructure team is not finished and it is designed to work with individuals names that are in the format POPINDNUMBER (eg France1, France2, Yoruba1)
		
		../../fs_4.0.1/fs_linux_glibc2.3 finestructure  -x 100000 -y 100000 -z 1000 output.chunkcounts.out structure_result.xml
		../../fs_4.0.1/fs_linux_glibc2.3 finestructure  -x 100000 -k 2 -m T -t 1000000 output.chunkcounts.out structure_result.xml structure_tree.out
		../../fs_4.0.1/fs_linux_glibc2.3 finestructure  -e meancoincidence  output.chunkcounts.out structure_result.xml structure_meancoincidence.csv
		../../fs_4.0.1/fs_linux_glibc2.3 finestructure  -e X2  output.chunkcounts.out structure_result.xml structure_meanstate.csv

9.-  To run GLOBETROTTER we need to run again chromopainter but this time not ALLvsALL but only with the query target population (eg SoAm-3) vs the putative donor populations (eg. PB1-Os PB1-ARG PA PB2 PB3 etc). Ideally the populations should be based now on the population clusters obtained from finestructure

10.-	Run chromopainter (for each chromosome) with an idfile with new population clusters (-t) and a file with target/donor populations (-f)
		
		for i in {1..16}; do ./ChromoPainterv2 -g chr${i}.haps -r chr${i}_rec -t idfile.txt -f CL210/CL210_poplist.txt 0 0 -o CL210/chr${i}; done &

11.-	Run GLOBETROTTER with chromopainter results. Create a parameter file, a samplelist file (text with the location of the files generated by chromopainter), and a recombination file list (text with the location of the recombination files generated in step 5)

		R < GLOBETROTTER.R PA/test_PA.params.txt PA_samplelist.txt reclist.infile.txt --no-save > testPA.txt &

		Run first with the option in parameter file null.ind: 1 (this generates a result indicating if there is detectable admixture or not). This wont work if the target population has only one individual!
		In case of admixture, run with null.ind : 0
		Perform boostrapping (only if you have more than one individual in target population)



INSTALLATION OF GSL version 1.6 for FINESTRUCTURE
sudo apt-get install build-essential
tar -xzvf /mnt/c/ubuntu/seubayanus_resubmission/beagle/../../fs_4.0.1/fs_linux_glibc2.3_4.0.1/gsl-1.6.tar.gz 
cd gsl-1.6/
./configure --prefix=$HOME/bin
make
make install
export GSL_CONFIG=$HOME/bin/bin/gsl-config
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/bin/lib
