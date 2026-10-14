# Day 3 Practical - 01. K-mer analyses

This first part of the practical will cover the steps for K-mer analyses for sex chromosome discovery. We will largely follow the **[KmersGwas](https://github.com/voichek/kmersGWAS/tree/master)** pipeline.  
This builds on the program **[KMC](https://github.com/refresh-bio/KMC)** to generate Kmer tables and **[PLINK](https://www.cog-genomics.org/plink/2.0/assoc)** for assocation tests.  

Start FileZilla with today's IP, open your terminal and connect to the server using your key file and username.  
Set up a new working directory for today, and activate the conda environment.

## 00. Prepare work folder for day 3

```
ssh -i ~/YOURLOCALFOLDER/user.pem user@44.244.12.10
mkdir day3
cd day3
# load conda environment and specify path to kmersgwas scripts
conda activate /opt/conda-envs/kmersgwas
export LD_LIBRARY_PATH=/opt/conda-envs/kmersgwas/lib:$LD_LIBRARY_PATH
# Software is located at:
/opt/course-software/kmersGWAS
```

## 01. Setup for K-mer analyses

Within day3, set up directories for K-mer analyses, copy all scripts needed and copy the phenotype data to your directory

```
mkdir kmersGWAS
cd kmersGWAS
mkdir Simmrg
cd Simmrg
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_phenotype.list .
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_MALE.list .
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_FEMALE.list .
```

## 02. Generate K-mer counts

This first step generates a list of all kmers and their presence/absence across all individuals. It can be run with the script below, but takes a very long time, so we won't run it today. Have a look at the file using the code below, and see if you can understand what it does.

```
cat ~/Share/day3/kmersGWAS/scripts/run_kmersGWAS_step1_Simmrg.sh
```

The output of this step is a file called *kmers.table. This is a long list of every kmer found in the samples and their presence/absence in each individual. The file is already saved in the server, and we can now use this file for the later steps.  This file can also be useful for other tools in the KMC/kmersGWAS tools, such as kinship analysis. Here we want to use this table to find kmers associated with sex, based on their presence or absence in our males and females.  

```
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_kmers_table* .
```

## 03. Test for association of K-mers and phenotype    
Testing for association with sex using the software **[PLINK](https://www.cog-genomics.org/plink/)**   
We used PLINK already for the GWAS yesterday, and here it is working in a very similar way.  
First generate PLINK compatible input file and filter for allele frequency:  
Use the script copied before, and set the input file, the kmer size, the phenotype file, minor allele frequency, minor allele count and output file name. Use the following parameters:  

t : input file (kmers_table)  
k : 31 (this is a standard k for species barcoding, and works well for delimiting sexes as well). This was already set in step1, but here we need to tell PLINK what we used.  
p : Sex phenotype list copied before    
maf : minor allele frequency, set to 0.05, standard filtering  
mac : minor allele count. Another filtering set, based on the number of times the minor allele occurs, rather than frequency. In a large dataset, this can be set to 5, as a lower limit, but if using a smaller dataset, where there are <5 samples per sex, this needs to be adjusted for the data. Here we have 5 males and 7 females, so we set mac to 4, to strictly filter.  
b : batch count. If the output is larger than b bytes, split into multiple files. Set to 1000000000 (1Gb). If multiple batches are created, each will end with .0, .1 etc.  
o : output name  

```
~/Share/day3/kmersGWAS/scripts/kmersGWAS/bin/kmers_table_to_bed -t kmers_table -k k -p sex_phenotypefile --maf maf --mac mac -b b -o output_name
```
If this doesn't work, you can copy the output: 
```
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_kmerGWAS_plink* .
```
Then run association test with PLINK. Set the following parameters:  
noweb : Don't check for online updates, like yesterday  
bfile : input file, output from the previous step (use your output name, plus the batch number ending '.0', without the ending '.bed')  
allow-no-sex : no sex as co-variate, like yesterday  
assoc : test for an association  
out : output file name  

```
plink --noweb --bfile plink_file --allow-no-sex --assoc --out outname
```

WARNING: this step may fail due to lack of memory. This shouldn't be an issue when working in your own university/institute cluster. If it does, copy the output from the shared folder.

```
cp ~/Share/day3/kmersGWAS/Simmrg/Simmrg_kmers.assoc .
```

This file is not tab-delimited, so annoying to work with. Use [awk](https://speciationgenomics.github.io/awk/) to fix this, and add the suffix .tab.  

```
awk -v OFS='\t' '{ $1=$1; print }' Simmrg_kmers.assoc > Simmrg_kmers.assoc.tab
```

This output file now lists every kmer and how associated it is with the phenotype, using a p-value.

You can download the output for this in R and see what p-value to set, but it is very large, and can take a long time. For simplicity, we will use only the smallest p-value, so only the most significantly associated kmers are filtered for (this is recommended by the authors of SexFindR).  

Check the output file for the most significant p-value, and filter for only kmers with this value.  
First cut out the 9th column (*cut -f 9*), then remove the header line containing 'P' (*grep -v 'P'*), then sort numerically (*sort -g*), and finally give the last value in this list (*head -1*).

```
cut -f 9 Simmrg_kmers.assoc.tab | grep -v 'P' | sort -g | head -1
```

This might take some time, and should results in the value 9.634e-07.  
Now use this value to filter the list of kmers, using awk.  
Edit the command for the p-value output from the step above.  

```
awk '$9 <= P' Simmrg_kmers.assoc.tab > Simmrg_most_significant_assoc.tab
```

Now, we have a list of kmers that are strongly associated with sex, but no information about where these kmers came from. The next step is to figure this out, using blast.

## 04. Assemble contigs from significant kmers    

This step makes small contigs from the associated kmers we just found. This is useful for two reasons:  
1. We are looking for a region of the genome associated with sex, so if there are many kmers associated with sex coming from the same region, this is more likely to be a sex-determining region. If many kmers align together and assemble into single contigs, this is highly informative of a large area of many associated kmers, therefore a good candidate region.
2. Our kmers are small, 31bp, and there are more likely to be misalignments of shorter nucleotide sequences in blast. Assembling into larger contigs makes the blast more reliable.

This uses a custom script from the kmerGWAS tool to convert the PLINK association file to a fasta file of all kmers. This can be used for assembly of potentially sex-specific contigs using the short read assembler **[ABySS](https://github.com/bcgsc/abyss)**

```
python ~/Share/day3/kmersGWAS/scripts/plink_to_abyss_kmers.py Simmrg_most_significant_assoc.tab Simmrg_abyss_input.txt
```

Assemble small contigs with ABySS. Use the following arguments:  
k : the span of a k-mer pair for ABySS. Shorter that the kmer size, so the kmer pairs are *mostly* overlapping. Set to 25.  
c : Coverage threshold. Remove contigs with mean k-mer coverage less than this value. Set to 0.  
e : Erosion threshold. Erode bases at the ends of blunt contigs with coverage less than this value. Set to 0.  
o : output name
```
ABYSS -k k -c c -e e Simmrg_abyss_input -o output_name
```

## 05. Locate the assembled contigs in reference genome    

Use **[BlastN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch)**  to locate contigs in the reference genome, for this first generate a blast reference database for the reference genome and then run blast.  

```
# make the reference database. Set the database type to nucleotide
makeblastdb -in ~/Share/day3/kmersGWAS/astcal.fasta -dbtype nucl
# run blast. Set the output format to #6, see the link for information.
blastn -query Simmrg_abyss_output.txt -db ~/Share/day3/kmersGWAS/astcal.fasta -outfmt 6 -out Simmrg_astcal_kmers_blast.out
```

## 06. Visualize genomic locations of the assembled contigs in R

Download the output file of the BLAST above Simmrg_astcal_kmers_blast.out file to your local machine. Also download the file chrlength_astcal.tsv from this ~/Share/day3/kmersGWAS .  
On your machine open **[R](https://cran.r-project.org)** or **[RStudio](https://posit.co/products/open-source/rstudio/?sid=1)**   

In R execute the code below. Note the blast filtering step, play around here and decide if you want to filter the output.  

```
install.packages("tidyverse")
library(tidyverse)

# Read in data
kmerblast_out <- read_tsv("Simmrg_astcal_kmers_blast.out") %>%
  rename("query_sid" = 1,
         "scaf" = 2,
         "ident.match_p" = 3,
         "alignmt_len" = 4,
         "mismatch_n" = 5,
         "gapopen_n" = 6,
         "query_start" = 7,
         "query_end" = 8,
         "ref_start" = 9,
         "ref_end" = 10,
         "exp_value" = 11,
         "bitscore" = 12)

# Read in reference chromosome list

ref <- read_tsv("chrlength_astcal.txt") %>%
  rename("scaf" = 1,
         "length" = 3) %>%
  filter(grepl('LS', scaf))
  

# have a look at the different variables in the blast output. Maybe you want to filter the blast hits based on this?

ggplot(kmerblast_out) +
  geom_histogram(aes(x = alignmt_len, y = after_stat(density)), binwidth = 1)
# + xlim(0, 200)

ggplot(kmerblast_out) +
  geom_histogram(aes(x = gapopen_n, y = after_stat(density)), binwidth = 1)

# query' is unique for each k-mer, so we want to select the best hit for each one and plot those:
kmerblast_filter <- 
  filter(kmerblast_out, grepl('LS', scaf)) %>%
# add extra filtering based on the histograms above
  group_by(query_sid) %>% 
  top_n(n = 1, wt = -exp_value)

# Count hits per chromosome

kmerblast_chrom <-
  kmerblast_filter %>% 
  ungroup() %>% 
  count(scaf)

kmer_chr <- left_join(ref, kmerblast_chrom) %>% arrange(scaf)

kmer_chr$n_prop <- kmer_chr$n / kmer_chr$length

(chr <-
    ggplot(kmer_chr,
           aes(x=scaf,
               y=n_prop)) + 
    geom_col() +
    theme_bw() +
    labs(title = bquote("Blasted kmerGWAS results for P. marginatus (North)"),
         subtitle = "Counts per chromosome, proportional to chromosome length") +
    xlab("Chromosome") +
    ylab("ABySS contig blastn hit proportion") +
    theme_bw() + 
    theme(legend.position = "none") + 
    theme(axis.text.x = element_text(size = 15, vjust = .5,  angle = 90),
          axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 10),
          axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(size=15)
    )
)
  
### kmer positions
  
  (points <-
      filter(kmerblast_filter, scaf == "my_chrom_name") %>%
      ggplot(aes(x = ref_start, y = ident.match_p)) +
      geom_point(
        size = 1) +
      labs(title = bquote("Blasted kmerGWAS results for P. marginatus (North)"),
           subtitle = "My chromosome")) +
      labs(x = "Position (bp)", y = "Blast similarity best hit per contig (%)") +
      theme_bw() + 
      theme(legend.position = "none") + 
      theme(axis.text.x = element_text(size = 15, vjust = .5,  angle = 90),
            axis.title.x = element_text(size = 15),
            axis.text.y = element_text(size = 10),
            axis.title.y = element_text(size = 12, angle = 90, hjust = .5, vjust = .5, face = "plain"),
            title = element_text(size=15)
      )

```
