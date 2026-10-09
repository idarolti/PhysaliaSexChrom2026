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
```

## 02. Run SCINKD3 workflow
Scinkd3 is a workflow that uses **[SnakeMake](https://snakemake.readthedocs.io/en/stable/)** , a workflow management system that can be used to generate analysis workflows via a human readable, Python based language. You can run the entire workflow by setting all needed information in a configuration file that is in **[json](https://www.json.org/json-en.html)** format.

### A. Understand the config file

Copy the config file to your working directory and inspect it.

```
mkdir SCINKD3
cd SCINKD3
cp /opt/course-software/SCINKD3/config_SCINKD.v3.1.5.json .
cat config_SCINKD.v3.1.5.json
```
You will see information towards the settings of the tool and the input data. This tool uses a reference genome in fasta format and several male and female short read data files in fastq.gz format.
Reference genome line 	"genome": "GCF_049243985.1_Ppicta_v3.0_genomic.fna")
This assembly has 18 chromosomes (line 	"ChrNum": "18")
Read file naming suffix for paired-end data 	"R1_suffix": "_1.fastq.gz", "R2_suffix": "_2.fastq.gz",

### B. Prepare the input data

We will run the tool for _Sphaerodactylus_ _townsendi_, Townsend's dwarf sphaero or Townsend's least gecko, and we will identify its sex chromosomes using male and female RAD data that have been generated for this purpose before (**[Pinto et al 2022](https://doi.org/10.1093/jhered/esac016)**).

The sequencing data from Pinto et al. 2022 are deposited under BioProject ID PRJNA746057
Inspect the BioProject **[here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA746057/)** to familiarise yourself with such data.

Then, download the reference genome from NCBI for _Sphaerodactylus townsendi_ (GCF_021028975.2) and unzip it
For this, you can search the NCBI Genome database, then chose the FTP site and copy the full path. In your server session download the file and unzip. Place it in the day3/scinkd3 directory
This assembly has 17 chromosomes

<details>

<summary>Prepare Genome</summary>

```
cd ~/day3/scinkd3
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/028/975/GCF_021028975.2_MPM_Stown_v2.3/GCF_021028975.2_MPM_Stown_v2.3_genomic.fna.gz
gunzip GCF_021028975.2_MPM_Stown_v2.3_genomic.fna.gz
```

</details>

Get the male and female short read data, we will use a subset of three individuals to run the analysis.
You could normally fetch those data from the BioProject mentioned above with the sratoolkit using the commands prefetch + the SRR accession number followed by the command fastq_dump.
For the sake of time, we have prepared these already so you can simply copy them to your directory.

```
cp ~/Share/day3/scinkd3/SRR15111* .
```
The individuals SRR15111480 to SRR15111486 are males, the individuals SRR15111473 to SRR15111479 are females

Now modify the config file in the SCINKD3 directory to match the genome name, the prefix and the sample names. Change all memory settings to "2"

<details>

<summary>Solution Modify Config File</summary>
You can either download the file to your machine and modify it with a text editor of your choice or use nano in the terminal to change the text to the one below

```
nano SCINKD3/config_SCINKD.v3.1.5.json
{
	"repair_mem": 2,
	"meryl_mem": 2,

	"per_job_threads": 2,
	"depth_window": "50000",

	"genome": "GCF_021028975.2_MPM_Stown_v2.3_genomic.fna",
	"ChrNum": "17",

	"R1_suffix": "_1.fastq.gz",
	"R2_suffix": "_2.fastq.gz",

	"males": [
	"SRR15111480",
	"SRR15111481",
	"SRR15111482",
	"SRR15111483",
	"SRR15111484",
	"SRR15111485",
	"SRR15111486"],

	"females": [
	"SRR15111473",
	"SRR15111474",
	"SRR15111475",
	"SRR15111476",
	"SRR15111477",
	"SRR15111479"]

}

Ctrl+O
Ctrl+X
```

</details>



### C. Start the Run

With SnakeMake you can first launch a dry run which checks if all is in place to run the workflow
```
cd ~/day3/scinkd3
snakemake --snakefile /opt/course-software/SCINKD3/SCINKD.v3.1.5.snakefile --cores 2 -np 
```

If this looks fine, you could start the workflow with this
```
snakemake --snakefile /opt/course-software/SCINKD3/SCINKD.v3.1.5.snakefile --cores 2 
```
**CAUTION: Since there is only 2 cores per student, this will FAIL.**


Successfully started, snakemake will let you know what it is doing, the workflow would run for quite some time so in the sake of time and due to limited resources we will not run it.
SCINKD3 will run through the steps illustrated below, you can also see them with inspecting the snakefile, let's do this

```
cat /opt/course-software/SCINKD3/SCINKD.v3.1.5.snakefile
```
SCINKD3 builds on the Kmer tool **[Meryl](https://github.com/marbl/meryl)** to count kmers first in each sample followed by intersecting kmers between individuals of the same sex and comparing the combined female and male Kmer catalog to identify kmers specific to each sex. It then uses the Meryl function "meryl-lookup" to identify reads that contain sex-specific kmers and then places those reads onto the reference genome using the read aligner **[Minimap2](https://github.com/lh3/minimap2)**. Ultimately, the results are presented as read coverage per sex of sex-specific reads along the genome as well as a plot that depicts the chromosomes with an accumulation of sex-specific reads. 

<img width="2070" height="1499" alt="Workflow" src="https://github.com/user-attachments/assets/72877682-4592-4bcc-af20-9456b50fdf4d" />

Meryl is quite memory intense and needs more resources then you have access to, so you can download the data we have run outside of the Amazon server.
SnakeMake has the functionality to pick up a run at intermediate steps, so this is what we will do here

```
cd ~/day3/scinkd3 
cp -r ~/Share/day3/scinkd3/SCINKD3RESULTS/* .
snakemake --snakefile /opt/course-software/SCINKD3/SCINKD.v3.1.5.snakefile --cores 2 
```
This should generate the output plots. If this is not working for you, you can also find the plots already prepared and download them in the next step.

### D. Inspect the output
The SnakeMake run should have produced the final plots in your day3/sckind3 directory. Else, uou can find the files SCINKD3.males.png, SCINKD3.females.png and SCINKD3.dotplot.png here
```
~/Share/day3/scinkd3/outputplots
```
Download these files and inspect them.
Which chromosome looks like a sex chromosome and what type of heterogamety do we have here?

<details>

<summary>Final Plots</summary>

<img width="1600" height="600" alt="SCINKD3 males" src="https://github.com/user-attachments/assets/b7068d7d-5eb7-40d9-9a8f-ea74a5cd8f12" />
<img width="1600" height="600" alt="SCINKD3 females" src="https://github.com/user-attachments/assets/3fc36fd2-22a3-4957-8954-b71a8698c022" />

<img width="1800" height="750" alt="SCINKD3 dotplot" src="https://github.com/user-attachments/assets/5b66663f-5efb-4745-9205-3158e06aad8d" />

You can search NCBI for the chromosome name "NC_059427.1" and you will see that this is chromosome 3. If you check in the original paper (**[Pinto et al 2022](https://doi.org/10.1093/jhered/esac016)**), you will see that they also had found chromosome 3 as a XY chromosome system.
</details>







