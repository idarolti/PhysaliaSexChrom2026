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
We will explore how this tool works across species with different degrees of sex chromosme divergence
Select to work on one of the following species.  
**[Vulpes vulpes](https://www.inaturalist.org/taxa/42069-Vulpes-vulpes)** - Red Fox.  
**[Anniella stebbinsi](https://www.inaturalist.org/taxa/479459-Anniella-stebbinsi)** - Southern California or San Diegan legless lizard.  
**[Lepidodactylus listeri](https://www.inaturalist.org/taxa/34352-Lepidodactylus-listeri)** - Christmas Island gecko.  
**[Sphaerodactylus notatus](https://www.inaturalist.org/taxa/33689-Sphaerodactylus-notatus)** - Florida Reef Gecko.  
  
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

bgzip mVulVul1.hap1.fasta.gz
bgzip mVulVul1.hap2.fasta.gz

  ```
</details>


<details>
  <summary>Anniella stebbinsi</summary>


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

The data for the Geckos are not available from the INSDC but over FigShare from the SCINKD2 paper.
For simplicity we placed them for you in the folder for day5

<details>
  <summary>Lepidodactylus listeri</summary>

  ```
cp ~/Share/day5/SCINKD/genomes/LepLis/Lepidodactylus_listeri_hic.hap1.fasta .
bgzip LepLis.hap1.fasta
cp ~/Share/day5/SCINKD/genomes/LepLis/Lepidodactylus_listeri_hic.hap2.fasta .
bgzip LepLisic.hap2.fasta
  ```
</details>

<details>
  <summary>Sphaerodactylus notatus</summary>

  ```
cp ~/Share/day5/SCINKD/genomes/SphNot/S_notatus_TG4245_Omni-C_hap1.map.pretext.fasta .
mv S_notatus_TG4245_Omni-C_hap1.map.pretext.fasta SphNot.hap1.fasta
bgzip SphNot.hap1.fasta
cp ~/Share/day5/SCINKD/genomes/SphNot/S_notatus_TG4245_Omni-C_hap2.map.pretext.fasta .
bgzip Lepidodactylus_listeri_hic.hap2.fasta
  ```
</details>
