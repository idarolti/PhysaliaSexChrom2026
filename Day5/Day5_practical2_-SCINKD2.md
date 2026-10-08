# Day 5 Practical 2 - Sex Chromosome Identification by Negating Kmer Densities (SCINKD) of haplotype-phased assemblies

This practical will use the SCINKD2 workflow to analyse different genome assemblies for their presence of sex chromosomes.
The original code is here **https://github.com/DrPintoThe2nd/SCINKD**
The corresponding publication is **[Pinto et al. 2026 MBE](https://doi.org/10.1093/molbev/msag067)**


## 01. Prepare work folder for day 5 practical 2

Open your terminal application as usual, connect to the server as usual, replace "25" with the number of your user account and use today's IP address

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@35.89.239.29
cd day5
conda activate /opt/conda-envs/day5
```

Don't forget to open a Filezilla connection, and change the IP to today's address


## 02. Download genome assemblies
We will explore how this tool works across species with different degrees of sex chromosme divergence
Select to work on one of the following species. 
**[Anniella stebbinsi](https://www.inaturalist.org/taxa/479459-Anniella-stebbinsi)** - Southern California or San Diegan legless lizard. 
**[Vulpes vulpes](https://www.inaturalist.org/taxa/42069-Vulpes-vulpes)** - Red Fox. 
**[Lepidodactylus listeri](https://www.inaturalist.org/taxa/34352-Lepidodactylus-listeri)** - Christmas Island gecko. 
**[Sphaerodactylus notatus](https://www.inaturalist.org/taxa/33689-Sphaerodactylus-notatus)** - Florida Reef Gecko. 
  


Open a web browser, search **[NCBI Genomes](https://www.ncbi.nlm.nih.gov/genome/)** for the entry "Apeltes quadracus"  
Select assembly **[GCA_048569185.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_048569185.1/)**  
Download the **[X](https://www.ncbi.nlm.nih.gov/nuccore/CM109091.1?report=fasta)** and the **[Y](https://www.ncbi.nlm.nih.gov/nuccore/CM109092.1?report=fasta)** chromosome fasta files to your server account  
You can either download them with the links above to your local machine and then use FileZille or from within the server run curl to download the files and then change their filename to make them easier to recognise.

```
 curl -o CM109091.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM109091.1&rettype=fasta&retmode=text"
mv CM109091.1.fasta XChrom.fasta
 curl -o CM109092.1.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM109092.1&rettype=fasta&retmode=text"
mv CM109092.1.fasta YChrom.fasta

```
         
## 04a. Option 1 Align X and Y chromosome online
Align the X and Y chromosome sequence online with **[DGenies](https://dgenies.toulouse.inra.fr)**  

```
https://www.ncbi.nlm.nih.gov/genome/](https://dgenies.toulouse.inra.fr
```
Click on the "RUN" tab  


<img width="1293" height="831" alt="Screenshot 2025-09-29 at 17 12 08" src="https://github.com/user-attachments/assets/bbc23361-b8be-4c76-b17b-817a79c2941d" />

  
Select the fasta files you downloaded previously and fill in the required information.  
