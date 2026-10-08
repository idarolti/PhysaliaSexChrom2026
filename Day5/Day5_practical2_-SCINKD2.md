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
conda activate /opt/conda-envs/day5
```

Don't forget to open a Filezilla connection, and change the IP to today's address


## 02. Download and prepare genome assemblies
We will explore how this tool works across species with different degrees of sex chromosme divergence
Select to work on one of the following species.  
**[Anniella stebbinsi](https://www.inaturalist.org/taxa/479459-Anniella-stebbinsi)** - Southern California or San Diegan legless lizard.  
**[Vulpes vulpes](https://www.inaturalist.org/taxa/42069-Vulpes-vulpes)** - Red Fox.  
**[Lepidodactylus listeri](https://www.inaturalist.org/taxa/34352-Lepidodactylus-listeri)** - Christmas Island gecko.  
**[Sphaerodactylus notatus](https://www.inaturalist.org/taxa/33689-Sphaerodactylus-notatus)** - Florida Reef Gecko.  
  
The first two species have haplotype phased assemblies in the INSDC record, so you can search for them as in practical 1 of today, identify both haplotypes, copy the download link from the FTP page and download the files.

<details>
  <summary>Anniella stebbinsi</summary>
  
  ### Anniella stebbinsi
  ```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/312/545/GCA_051312545.2_rAnnSte1.2_hap2/GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/051/312/545/GCA_051312545.2_rAnnSte1.2_hap2/GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz

mv GCA_051312545.2_rAnnSte1.2_hap2_genomic.fna.gz rAnnSte1.2_.hap1.fasta.gz
mv GCA_051312515.2_rAnnSte1.2_hap1_genomic.fna.gz rAnnSte1.2_.hap2.fasta.gz

gunzip rAnnSte1.2_.hap1.fasta.gz
gunzip rAnnSte1.2_.hap2.fasta.gz

bgzip rAnnSte1.2_.hap1.fasta
bgzip rAnnSte1.2_.hap2.fasta
  ```
</details>


```
 curl -o CM109091.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM109091.1&rettype=fasta&retmode=text"
mv CM109091.1.fasta XChrom.fasta
 curl -o CM109092.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM109092.1&rettype=fasta&retmode=text"
mv CM109092.1.fasta YChrom.fasta

```
For easier handling we will rename the files, SCINKD2 also has the following restriction so we take care of this as well
File naming restriction: Both input haplotype fasta files MUST be bgzipped and MUST end in ".hap1.fasta.gz" and ".hap2.fasta.gz"

```
```


## 04a. Option 1 Align X and Y chromosome online
Align the X and Y chromosome sequence online with **[DGenies](https://dgenies.toulouse.inra.fr)**  

```
https://www.ncbi.nlm.nih.gov/genome/](https://dgenies.toulouse.inra.fr
```
Click on the "RUN" tab  


<img width="1293" height="831" alt="Screenshot 2025-09-29 at 17 12 08" src="https://github.com/user-attachments/assets/bbc23361-b8be-4c76-b17b-817a79c2941d" />

  
Select the fasta files you downloaded previously and fill in the required information.  
