# Day 2 Practical Sex chromosome discovery - coverage and SNP based methods

This practical will cover methods on sex chromosome discovery:

1. Coverage based analysis
2. SNP-based analyses
    
We will base our analysis on the **[SexFindR pipeline](https://sexfindr.readthedocs.io/en/latest/#)** published by **[Grayson et al](https://doi.org/10.1101/2022.02.21.481346)**


## 00. Prepare work folder for day 2
Open your terminal application as yesterday, connect to the server as yesterday, using your .pem file and user account, and use today's IP adress

```
ssh -i ~/YOURLOCALFOLDER/user.pem user@44.251.209.2
mkdir day2
cd day2
conda activate /opt/conda-envs/day2
```
Don't forget to open a Filezilla connection, and change the IP to today's address  

## 01. Coverage based analysis (Guppies)
We will follow the workflow of **[DifCover](https://github.com/timnat/DifCover)** to generate a coverage file for females and males, we start with the coordinate sorted BAM files generated yesterday.       
Set up directories and files

```
mkdir coverage
cd coverage
mkdir picta
```

### Step 1. Create BED file

Using **[BEDtools](https://bedtools.readthedocs.io/en/latest/)** we transform the coverage information of the BAM files into a unionbed file   
This step takes a long time, so start running to check it works, then cancel (ctrl+c).

```
cd picta
~/Share/day2/coverage/dif_cover_scripts/from_bams_to_unionbed.sh ~/Share/day2/bam_files/Poecilia_picta_female1_subset.bam ~/Share/day2/bam_files/Poecilia_picta_male2_subset.bam
```

### Step 2. Calculate the coverage ratio per window

In order to run this command, use the premade output of day2.
We will calculate average coverage of valid bases across all merged bed intervals for the female and the male file

Inspect the coverage statistics from the covstats file located in Day2/coverage to determine the parameter values for P. picta for the command below for P. picta.

```
cat ~/Share/day2/coverage/covstats.tab
```

Use these values to set the parameters for the following command

a = minimum coverage sample1  
A = maximum coverage sample1  
b = minimum coverage sample2  
B = maximum coverage sample2  
v = target number of valid bases in the window (set to 10000)  
l = minimum size of window to output (set to 1000)   

Then run the command below to generate the unionbed coverage file. This takes a few minutes.

```
~/Share/day2/coverage/dif_cover_scripts/from_unionbed_to_ratio_per_window_CC0 -a a -A A -b b -B B -v v -l l ~/Share/day2/coverage/picta_sample1_sample2.unionbedcv
```

### Step 3. Adjust coverage based on the bam ratio in covstats.tab and generate DNAcopy output file
The covstats.tab file also contains a value 'bam_ratio', which is the ratio of the sizes of the two .bam files. We use this to adjust the ratio per window, to account for the fact that one sample might have been sequenced to a higher coverage than the other. This then adjusts the ratio calculated in the previous step by this value. Change the value of 'bam_ratio' in the following command to the value in covstats.tab, and run.  
If step 2 took to long, retrieve the input for step3 like this, if not proceed to second command
```
cp ~/Share/day2/coverage/picta_sample1_sample2.ratio_per_w_CC0_a0_A20015_b0_B169426_v10000_l1000 ./sample1_sample2.ratio_per_w_CC0_a0_A20015_b0_B169426_v10000_l1000
```

```
~/Share/day2/coverage/dif_cover_scripts/from_ratio_per_window_to_prepare_for_DNAcopy_output.sh sample1_sample2.ratio_per_w_CC0_* bam_ratio
```

### Step 4. Create coverage plots in R

This script uses the R package **[DNAcopy[(https://bioconductor.org/packages/release/bioc/html/DNAcopy.html)** to generate output plots of the coverage ratio. DNAcopy uses both the raw coverage ratio, and computes the average coverage ratio of adjacent windows to create the plot. Run this on the server
```
Rscript ~/Share/day2/coverage/dif_cover_scripts/run_DNAcopy_from_bash.R sample1_sample2.ratio_per_w_CC0_*.log2adj_*
```
Dowload the pdf file to your machine using Filezilla (or equivalent) and inspect it.

### Step 5. Filter only genomic regions with enrichment scores > p.

The script extracts from file *.DNAcopyout fragments with enrichment scores ≥ p and stores them in *.DNAcopyout.up{p}, (i.e. fragments where read coverage in sample1 is higher than sample2 ), and *.DNAcopyout.down{-p} fragments with enrichment scores ≤-p, (i.e. fragments where coverage in sample2 is higher than sample1 ).

Set p to 2, to filter sites with double coverage in one sample compared to the other.

```
~/Share/day2/coverage/dif_cover_scripts/from_DNAcopyout_to_p_fragments.sh sample1_sample2.*.DNAcopyout p
```

#### Now run this again for P. reticulata  
#### If you finish both species, try adjusting the window size in step 2, and see how this changes the results. Try smaller windows, e.g. 1,000 or larger, e.g. 50,000


   
## 02. SNP based analyses (Cichlids)

If you logged out, remember to reactivate the conda environment

```
conda activate /opt/conda-envs/day2
```

On the server in your home set up directories and files.

We will use samples from the Lake Tanganyika cichlid Pseudosimochromis marginatus (North variant), Simmrg for short.

```
cd 
cd day2
mkdir SNPbased
cd SNPbased
mkdir Fst
mkdir GWAS
mkdir SNPden
cp ~/Share/day2/SNPbased/Simmrg_FEMALE.list Fst/
cp ~/Share/day2/SNPbased/Simmrg_MALE.list Fst/
cp ~/Share/day2/SNPbased/Simmrg_sex.list GWAS/
cp ~/Share/day2/SNPbased/Simmrg_FEMALE.list SNPden/
cp ~/Share/day2/SNPbased/Simmrg_MALE.list SNPden/
cp ~/Share/day2/SNPbased/Simmrg_astcal.final.vcf.gz .
```

**NOTE** If any steps fail, the relevant files can be copied from ~/Share/day2/SNPbased

### Analysis 1 Calculate intersex Fst 
We will use **[VCFtools](https://vcftools.github.io)** to calculate Fst (fixation index, a measure of genetic differentiation between populations) between males and female in 10kb windows based on the variant file you learned how to generate yesterday   

Use the command vcftools with the following arguments:

gzvcf : the input VCF file (Simmrg_astcal.final.vcf.gz)    
weir-fst-pop : a list of samples in each population. Use this command twice for two populations (i.e. a list of male and female samples)    
fst-window-size : size of windows to calculate Fst in (use 10000)    
out : name of output file (i.e. Simmrg)    

```
cd Fst
vcftools --gzvcf vcf --weir-fst-pop sample_list1 --weir-fst-pop sample_list2 --fst-window-size windowsize --out outputname
```

#### Plot in R  
Download the output file and the file astcal.fasta.fai. This in an index of the reference this data was mapped to, and contains the lengths of each chromosome.
Download the Simmrg_Fst_plots.R script from day2/SNPbased and run in R.
Remember to set your working directory to where the files are saved.  

### Analysis 2 Run association test for sex with GEMMA   

We will first generate a filtered input file that we will then further format with PLINK and then use as input for **[GEMMA](https://github.com/genetics-statistics/GEMMA)** which runs linear mixed models to test for an association between each SNP and sex.

This starts again with vcftools to create a plink file based on filtering of the VCF file. Run with the following arguments:

gzvcf : the input VCF file (Simmrg_astcal.final.vcf.gz)       
plink : tells VCFtools to output a plink file    
remove-indels : keep only SNPs, not indels (insertions/deletions up to 50bp)    
max-missing : Set the maximum allowed number of samples with missing data at a site. This is not important in our case because the data has already been filtered to remove sites with any missing data, but good practice to set to 0.5 - we need sites present in both males and females.    
max-maf / maf : Include only sites with a Minor Allele Frequency greater than or equal to the "--maf" value and less than or equal to the "--max-maf" value. One of these options may be used without the other. Allele frequency is defined as the number of times an allele appears over all individuals at that site, divided by the total number of non-missing alleles at that site. Set --maf to 0.05 and --max-maf to 0.95.    
out : name of output file (i.e. Simmrg_step1)    

```
cd GWAS
vcftools --gzvcf vcf --plink --remove-indels --max-missing max-missing --max-maf max-maf --maf maf --out outname
```

Now we run plink to summarise the SNPs into values for each sex. This used the sex.list which lists samples in the 2 populations, males and females, with values 1 and 2. Use 'cat' to inspect this file.

Run with the following arguments:

file : your output from the previous step. Note there are 2 output files, but you can just write the outname used in the previous step, plink will find the file needed    
pheno : the list of samples, allocated to sexes    
make-bed : create a binary file set of the output with the filtering parameters stored. Good for logging.    
out : name of output file for this step    
no-web : don't check for plink updates    
allow-no-sex : We are using sex as a phenotype, and not a covariable in our data. Therefore we do not provide 'sex' information, and need to tell plink not to remove samples without sex information.     

```
plink --file vcf_output --pheno sex.list --make-bed --out outname --noweb --allow-no-sex
```

GEMMA performs the actual GWAS. Use the following commands:

bfile : output from plink step (without suffix)    
lm : type of linear model to use (1: Wald test, 2: Likelihood ratio test, 3: Score test, 4: all). Use Likelihood ratio test.    
o : output file name.    

```
gemma -bfile plink_output -lm lm -o outname
```

Result is in the output/ folder. Download the assoc.txt file and the Simmrg_GWAS_plots.R script from day2/SNPbased and run in R. You will also need the file astcal.fasta.fai

### Analysis 3 Calculate male and female SNP density
We will format the SNP variant file with **[bcftools](https://samtools.github.io/bcftools/bcftools.html)** to subset female and male entries into separate files  
We will then calculate SNP density in 10kb windows with VCFtools.

First, unzip the VCF file using the command 'gunzip', and move it to the SNPden directory.

We will use bcftools to make a mini VCF with only one individual. This is done in a loop, first reading the list of female samples, and naming them all as females, and then the same for the list of male samples.

```
# open the for loop, and name each object "SAMPLE" based on a line in the file "Simmrg_FEMALE.list"
for SAMPLE in `cat Simmrg_FEMALE.list`;
  do
# use bcftools to extract the sample (-s) from the VCF and store it in its own VCF file (-o). -a trims alternative alleles that are not called as genotypes.
  bcftools view -a -s ${SAMPLE} -o FEMALE_${SAMPLE}.vcf Simmrg_astcal.final.vcf
# compress the single-sample VCF
  bgzip -c FEMALE_${SAMPLE}.vcf > FEMALE_${SAMPLE}.vcf.gz
# index the single-sample VCF
  bcftools index FEMALE_${SAMPLE}.vcf.gz
# calculate the number of SNPs per 10000 bp, and store this to an output file for this sample.
  vcftools --gzvcf FEMALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out FEMALE_${SAMPLE}
# close the for loop
done

for SAMPLE in `cat Simmrg_MALE.list`;
  do
  bcftools view -a -s ${SAMPLE} -o MALE_${SAMPLE}.vcf Simmrg_astcal.final.vcf
  bgzip -c MALE_${SAMPLE}.vcf > MALE_${SAMPLE}.vcf.gz
  bcftools index MALE_${SAMPLE}.vcf.gz
  vcftools --gzvcf MALE_${SAMPLE}.vcf.gz --SNPdensity 10000 --out MALE_${SAMPLE}
done
```

Copy all files called *.snpden to your local machine  
Open Rstudio, download the script Simmrg_SNPden_plots.R and load it in R. You will again need the reference index file.

#### Now do the same for Petrochromis polyodon (Petpol)
