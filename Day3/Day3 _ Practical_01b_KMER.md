# Day 3 Practical - 01b. K-mer analyses

This second part of the practical will use another tool building on K-mers analyses for sex chromosome discovery. 
We will use Sex Chromosome Identification by Negating Kmer Densities (SCINKD) for reference genomes and short reads.
 
Pinto BJ, Gable SM, Keating SE, Smith CH, Gamble T, Nielsen SV, Wilson MA. (2026). Sex chromosome identification and genome curation from a single individual with SCINKD. 
Molecular Biology and Evolution. 43(4). **https://doi.org/10.1093/molbev/msag067**

Start FileZilla with today's IP, open your terminal and connect to the server, set up a new working directory for today, and activate the scinkd3 conda environment

## 00. Prepare work folder for day 3

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@18.237.201.56
mkdir day3
cd day3

```

## 01. Setup for SCINKD3 analyses identifying sex chromosomes from short reads by comparing male and female catalogs 
We will first run scinkd3, which compares kmer catalogs of males and females using a SnakeMake Workflow

Within day3, set up a directory for scinkd3-based analyses and set the right compute environment 

```
mkdir scinkd3
cd scinkd3
conda activate /opt/conda-envs/day5
/opt/course-software/SCINKD3

```

## 02. Run SCINKD3 workflow
Scinkd3 is a workflow that uses **[SnakeMake](https://snakemake.readthedocs.io/en/stable/)** , a workflow management system is that can be used to generate analysis workflows via a human readable, Python based language. You can run the entire workflow by setting all needed information in a configuration file that is in **[json](https://www.json.org/json-en.html)** format.
Copy the config file to your working directory and inspect it.

```
cp ~/Share/day3/scinkd3/config_SCINKD.v3.1.5.json
cat config_SCINKD.v3.1.5.json
```
You will see information towards the settings of the tool and the input data.
Reference genome line 	"genome": "GCF_049243985.1_Ppicta_v3.0_genomic.fna")
This assembly has 18 chromosomes (line 	"ChrNum": "18")
Read file naming suffix for paired-end data 	"R1_suffix": "_1.fastq.gz", "R2_suffix": "_2.fastq.gz",


We will run the tool for another species

nd we will identify its sex chromosomes using male and female RAD data that have been generated for this purpose before ([**Pinto et al 2022](https://doi.org/10.1093/jhered/esac016)**)
Sphaerodactylus townsendi, Townsend's dwarf sphaero or Townsend's least gecko
this assembly has 17 chromosomes (line 	"ChrNum": "18")


Download the reference genome from NCBI for Sphaerodactylus townsendi and unzip it

```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/049/243/985/GCF_049243985.1_Ppicta_v3.0/GCF_049243985.1_Ppicta_v3.0_genomic.fna.gz
gunzip GCF_049243985.1_Ppicta_v3.0_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/028/975/GCF_021028975.2_MPM_Stown_v2.3/GCF_021028975.2_MPM_Stown_v2.3_genomic.fna.gz
gunzip GCF_021028975.2_MPM_Stown_v2.3_genomic.fna.gz
```
This first step generates a list of all kmers and their presence/absence across all individuals. It can be run with the script below, but takes a very long time, so we won't run it today. Have a look at the file using the code below, and see if you can understand what it does.

The RAD seq data from Pinto et al. 2022 are deposited under BioPorject ID PRJNA746057
Inspect the BioProject **[here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA746057/)** 

<details>

<summary>Tips for collapsed sections</summary>

### You can add a header

You can add text within a collapsed section.

You can add an image or a code block, too.

```ruby
   puts "Hello World"
```

</details>

The output of this step is a file called kmers.table.table. This is a long list of every kmer found in the samples and their presence/absence in each individual. The file is already saved in the server, and we can now use this file for the later steps.  

```
cp ~/Share/day3/kmersGWAS/picta/picta_kmers_table* .
```

## 03. Generate kinship table  

Generate kinship table, in case useful for future analysis:

```
../emma_kinship_kmers -t picta_kmers_table -k 31 --maf 0.05 > picta_kmers_table.kinship
cat picta_kmers_table.kinship
```

## 04. Test for association of K-mers and phenotype    
Testing for association with sex using the software **[PLINK](https://www.cog-genomics.org/plink/)**   
First generate PLINK compatible input file and filter for allele frequency

```
../kmers_table_to_bed -t picta_kmers_table -k 31 -p Ppicta_phenotype.txt --maf 0.05 --mac 2 -b 1000000000 -o picta_kmerGWAS_plink
```

Then run association test with PLINK

```
plink --noweb --bfile picta_kmerGWAS_plink.0 --allow-no-sex --assoc --out picta_kmers
```

WARNING: this step may fail due to lack of memory. This shouldn't be an issue when working in your own university/institute cluster. If it does, copy the output from the shared folder.

```
cp ~/Share/day3/kmersGWAS/picta/picta_kmers.assoc .
```

Edit format of output file

```
awk -v OFS='\t' '{ $1=$1; print }' picta_kmers.assoc > picta_kmers.assoc.tab
```

Check the output file for the most significant p-value, and filter for only kmers with this value

```
cut -f 9 picta_kmers.assoc.tab | grep -v 'P' | sort -g | head -1
```

This might take some time, and should results in the value 0.000532  
Edit the command for the p-value output from the step above.

```
awk '$9 <= P' picta_kmers.assoc.tab > picta_most_significant_assoc.tab
```

## 05. Assemble contigs from significant kmers    

Convert PLINK association file to a fasta input file for assembly of potentially sex-specific contigs using the short read assembler **[ABYSS](https://github.com/bcgsc/abyss)**

```
python ../plink_to_abyss_kmers.py picta_most_significant_assoc.tab picta_plink_abyss_input.txt
```

Assemble small contigs with ABYSS

```
ABYSS -k25 -c0 -e0 picta_plink_abyss_input.txt -o picta_plink_abyss_output.txt
```

## 06. Locate the assembled contigs in reference genome    

Use **[BlastN](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch)**  to locate contigs in the reference genome, for this first generate a blast reference database for the reference genome and then run blast

```
cp ~/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna .

makeblastdb -in Poecilia_picta.fna -dbtype nucl

blastn -query picta_plink_abyss_output.txt -db Poecilia_picta.fna -outfmt 6 -out picta_kmers_blast.out
```

## 07. Visualize genomic locations of the assembled contigs in R

Download the output file of the BLAST above picta_kmers_blast.out file to your local machine. Also download the file chrlength_Ppicta.tsv from this github.  
On your machine open **[R](https://cran.r-project.org)** or **[RStudio](https://posit.co/products/open-source/rstudio/?sid=1)**   

In R execute the code below

```
install.packages("tidyverse")
library(tidyverse)

# Read in data
kmerblast_out <- read_tsv("picta_kmers_blast.out") %>%
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

ref <- read_tsv("chrlength_Ppicta.tsv") %>%
  rename("scaf" = 1,
         "length" = 2)

# query' is unique for each k-mer, so we want to select the best hit for each one and plot those:
kmerblast_filter <- 
  filter(kmerblast_out, !grepl('JAVY', scaf)) %>% 
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
    labs(title = bquote("Blasted kmerGWAS results for P. picta"),
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
      filter(kmerblast_filter, scaf == "CM065364.1") %>%
      ggplot(aes(x = ref_start, y = ident.match_p)) +
      geom_point(
        size = 1) +
      labs(title = bquote("Blasted kmerGWAS results for P. picta"),
           subtitle = "All kmers")) +
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

