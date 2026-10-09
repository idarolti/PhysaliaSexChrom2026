# Day 5 Practical 2 - Sex Chromosome Identification by Negating Kmer Densities (SCINKD) of haplotype-phased assemblies

This practical will use the SCINKD2 workflow to analyse different genome assemblies for their presence of sex chromosomes.
The original code is here **https://github.com/DrPintoThe2nd/SCINKD**
The corresponding publication is **[Pinto et al. 2026 MBE](https://doi.org/10.1093/molbev/msag067)**
As on day3 this is a workflow that is executed over SnakeMake.

## 01. Prepare work folder for day 5 practical 2

Open your terminal application as usual, connect to the server as usual, replace "25" with the number of your user account and use today's IP address

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@35.89.239.29
cd day5
mkdir SCINKD2
cd SCINKD2
conda activate /opt/conda-envs/day5
```

Don't forget to open a Filezilla connection, and change the IP to today's address


## 02. Download and prepare genome assemblies
We will explore how this tool works across species with different degrees of sex chromosome divergence.  
Select to work on one of the following species.  
**[Vulpes vulpes](https://www.inaturalist.org/taxa/42069-Vulpes-vulpes)** - Red Fox (17 chromosomes).  
**[Anniella stebbinsi](https://www.inaturalist.org/taxa/479459-Anniella-stebbinsi)** - Southern California or San Diegan legless lizard (10 chromosomes).  
**[Lepidodactylus listeri](https://www.inaturalist.org/taxa/34352-Lepidodactylus-listeri)** - Christmas Island gecko (22 chromosomes).  
**[Sphaerodactylus notatus](https://www.inaturalist.org/taxa/33689-Sphaerodactylus-notatus)** - Florida Reef Gecko (17 chromosomes).  
  
The first two species have haplotype phased assemblies in the INSDC record, so you can search for them as in practical 1 of today, identify both haplotypes, copy the download link from the FTP page and download the files.
SCINKD2 has the following requirements to accept the genome files:
File naming restriction: Both input haplotype fasta files MUST be bgzipped and MUST end in ".hap1.fasta.gz" and ".hap2.fasta.gz"

<details>
  <summary>Vulpes vulpes</summary>
<img width="1510" height="1190" alt="Vulpesvulpes1" src="https://github.com/user-attachments/assets/4066c7c3-98c3-4583-bbc3-e1341a8632ab" />
<img width="857" height="795" alt="Vulpesvulpes2" src="https://github.com/user-attachments/assets/1ea67f86-fcd8-495d-9885-2e30a9abff98" />
<img width="1166" height="432" alt="Vulpesvulpes3" src="https://github.com/user-attachments/assets/f2c0e5d0-ae5c-4246-9f9e-de1b450f14e6" />

  ```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/925/GCA_964106925.2_mVulVul1.hap2.2/GCA_964106925.2_mVulVul1.hap2.2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/964/106/825/GCA_964106825.2_mVulVul1.hap1.2/GCA_964106825.2_mVulVul1.hap1.2_genomic.fna.gz

mv GCA_964106825.2_mVulVul1.hap1.2_genomic.fna.gz mVulVul1.hap1.fasta.gz
mv GCA_964106925.2_mVulVul1.hap2.2_genomic.fna.gz mVulVul1.hap2.fasta.gz

gunzip mVulVul1.hap1.fasta.gz
gunzip mVulVul1.hap2.fasta.gz

bgzip mVulVul1.hap1.fasta
bgzip mVulVul1.hap2.fasta

  ```
</details>


<details>
  <summary>Anniella stebbinsi</summary>


  ```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/312/545/GCA_051312545.2_rAnnSte1.2_hap2/GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/312/545/GCA_051312545.2_rAnnSte1.2_hap2/GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz

mv GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz rAnnSte1.2.hap1.fasta.gz
mv GCA_051312515.2_rAnnSte1.2_hap1_genomic.fna.gz rAnnSte1.2.hap2.fasta.gz

gunzip rAnnSte1.2.hap1.fasta.gz
gunzip rAnnSte1.2.hap2.fasta.gz

bgzip rAnnSte1.2.hap1.fasta
bgzip rAnnSte1.2.hap2.fasta
  ```
</details>

The data for the Geckos are not available from the INSDC but over FigShare from the SCINKD2 paper.
For simplicity we placed them for you in the folder for day5

<details>
  <summary>Lepidodactylus listeri</summary>

  ```
cp ~/Share/day5/SCINKD2/genomes/LepLis/Lepidodactylus_listeri_hic.hap1.fasta .
mv Lepidodactylus_listeri_hic.hap1.fasta LepLis.hap1.fasta
bgzip LepLis.hap1.fasta
cp ~/Share/day5/SCINKD2/genomes/LepLis/Lepidodactylus_listeri_hic.hap2.fasta .
mv Lepidodactylus_listeri_hic.hap2.fasta LepLis.hap2.fasta
bgzip LepLis.hap2.fasta
  ```
</details>

<details>
  <summary>Sphaerodactylus notatus</summary>

  ```
cp ~/Share/day5/SCINKD2/genomes/SphNot/S_notatus_TG4245_Omni-C_hap1.map.pretext.fasta .
mv S_notatus_TG4245_Omni-C_hap1.map.pretext.fasta SphNot.hap1.fasta
bgzip SphNot.hap1.fasta
cp ~/Share/day5/SCINKD2/genomes/SphNot/S_notatus_TG4245_Omni-C_hap2.map.pretext.fasta .
mv S_notatus_TG4245_Omni-C_hap2.map.pretext.fasta SphNot.hap2.fasta
bgzip SphNot.hap2.fasta

  ```
</details>

## 03. Run analysis
We will here show the example for **_Vulpes vulpes_**.  
For this you need to modify the config file to your genome of interest, similar to the exercise on day3; set the cores to 1 in the config and give snakemake 2
Please again only use two cores.
Number of chromosomes for _Vulpes vulpue_ is 17

```
mkdir SCINKD
cp /opt/course-software/SCINKD/config.json SCINKD/
nano SCINKD/config.json
{
	"per_job_threads": 1,
	"per_job_memory": 1,
	"ChrNum": 17,

	"prefix": "mVulVul1"
}
Ctrl+O
Ctrl+X
```

  
**THE CLUSTER HAS NOT THE RESOURCES TO RUN THIS SO INSTEAD DO THIS**
As on day3, first do a dry run.
We have prepared the intermediate steps for you, do the following to get those and to generate the final plots. Snakemake can restart from an interruppted workflow.
Try the dry run, then copying should take around 10 minutes, then redo the command as true command removing "-np"

```
time snakemake --use-conda --rerun-incomplete --nolock --cores 2 -s /opt/course-software/SCINKD/SCINKD.v2.2.4.snakefile -np #dry run
cp -r ~/Share/day5/SCINKD2/scinkd2_VulVul/RESULTS/* .
time snakemake --use-conda --rerun-incomplete --nolock --cores 2 -s /opt/course-software/SCINKD/SCINKD.v2.2.4.snakefile
```

For the next step we also need visualisation so we will once more generate pairwise alignments with Minimap2
```
minimap2 -x asm5 -t1 -c --eqx --secondary=no rAnnSte1.2.hap1.fasta.gz rAnnSte1.2.hap2.fasta.gz > rAnnSte1.2.asm.hic.paf
minimap2 -x asm5 -t1 -c --eqx --secondary=no mVulVul1.hap1.fasta.gz mVulVul1.hap2.fasta.gz > mVulVul1.asm.hic.paf
minimap2 -x asm5 -t1 -c --eqx --secondary=no LepLis.hap1.fasta.gz LepLis.hap2.fasta.gz > LepLis.asm.hic.paf
minimap2 -x asm5 -t1 -c --eqx --secondary=no SphNot.hap1.fasta.gz SphNot.hap2.fasta.gz > SphNot.asm.hic.paf
```

If the command above fails, you can copy those files from here to your local machine

```
~/Share/day5/SCINKD2/PAFFILES/ ##chose the paf file for your species
~/Share/day5/SCINKD2/
```

Can you identify the sex chromosomes from the plot?
Compare your plot files to those in **[Pinto et al. 2026 MBE](https://doi.org/10.1093/molbev/msag067)**


## 04. Plot data
Transfer the  files with the following file endings to your local machine: .results, .fai, .bed from here

```
~/Share/day5/SCINKD2/PAFFILES/ ##chose the paf file for your species
##chose your results files from one of the four species
~/Share/day5/SCINKD2/scinkd2_VulVul/RESULTS/
~/Share/day5/SCINKD2/scinkd2_AnnSte/RESULTS/
~/Share/day5/SCINKD2/scinkd2_LepLis/RESULTS/
~/Share/day5/SCINKD2/scinkd2_SphNot/RESULTS/
```


Open RStudio and load the script PlotSCINKD2.R from this repository


