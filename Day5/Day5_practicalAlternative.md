# Day 5 Practical

This practical will cover:

1. Generating whole reference genome alignments
2. Visualising dot plots
3. Investigate structural rearrangements

We will investigate sex chromosome inversions identified in the fourspine stickleback by **[Liu et al. 2025](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1011465)**


## 01. Local installation of visualisation software

* **[JBrowse](https://jbrowse.org/jb2/)** - The next-generation genome browser

Download JBrowse, a standalone application, to your local machine, chose the apropriate version for your system platform,
double click the installation and follow the on screen instructions, double click the JBrowse App Icon

```
https://jbrowse.org/jb2/download/
```


## 02. Prepare work folder for day 5

Open your terminal application as yesterday, connect to the server as yesterday, replace "25" with the number of your user account and use today's IP adress

```
ssh -i ~/YOURLOCALFOLDER/chrsex5.pem user5@44.251.209.2
mkdir day5
cd day5
conda activate /home/ubuntu/miniconda3/envs/sexchr
```

Don't forget to open a Filezilla connection, and change the IP to today's address


## 03. Identify sex chromosome assemblies from INSDC
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
Click **Submit**  

Your status will change from  

  
<img width="817" height="362" alt="Screenshot 2025-09-29 at 17 17 08" src="https://github.com/user-attachments/assets/3a5cab80-3197-4b2d-9ef4-7472a95e4dc1" />  


To  

<img width="785" height="460" alt="Screenshot 2025-09-29 at 17 18 53" src="https://github.com/user-attachments/assets/0432353d-228d-41eb-8bf1-07703a0b0e0f" />

in a couple of minutes<br/>  
<br/>

## 05a. Inspect the results online

<img width="787" height="746" alt="Screenshot 2025-09-30 at 09 51 35" src="https://github.com/user-attachments/assets/1b957ff2-a4a2-4428-b35c-3cd015cde012" />  
<br/>

## 06a. Export the results
Export the results in PAF format to your local machine.  
<br/>
<br/>

## 04b. Option 2 Align X and Y chromosome on the server
In your day 5 folder with the fasta files, align the two files with **[Minimap2](https://github.com/lh3/minimap2)**
```
minimap2 -c XChrom.fasta YChrom.fasta > XYAlign.paf
```
Transfer the file to your local machine using FileZilla and load it into Jbrowse

<br/>

## 07. Load the results into Jbrowse
Open Jbrowse  
Click **OPEN NEW GENOME**  

<img width="601" height="579" alt="Screenshot 2026-10-05 at 16 23 52" src="https://github.com/user-attachments/assets/f4ffc2d5-a70d-4658-80bf-0246a91919e8" />


Fill in the information, select your X chromosome fasta file and click **ADD ANOTHER ASSEMBLY**  at the bottom, fill in the information
  
<img width="603" height="632" alt="Screenshot 2026-10-05 at 16 24 42" src="https://github.com/user-attachments/assets/2f40f0cd-ee2d-4352-b33e-fadf7e004660" />

Click **SUBMIT**  
<br/>
Select  **Dotplot View -> LAUNCH VIEW**  
The following error message will appear  
<img width="1390" height="416" alt="Screenshot 2026-10-05 at 16 25 29" src="https://github.com/user-attachments/assets/1db1fe95-a027-4664-9572-5ee60d60556a" />
<br/>
Choose the following settings  
<img width="1399" height="574" alt="Screenshot 2026-10-05 at 16 25 51" src="https://github.com/user-attachments/assets/d98f31aa-7a32-4ec4-af26-6e65e704885d" />



Click **LAUNCH**  
and explore the plot
<br/>
Click **ADD** , choose **Linear Synteny View**  

<img width="1408" height="253" alt="Screenshot 2026-10-05 at 16 31 07" src="https://github.com/user-attachments/assets/9a08e2a8-da29-47b1-bb87-24924a57f385" />

Click **LAUNCH**  

You can now explore this plot  

<img width="1375" height="334" alt="Screenshot 2025-09-30 at 10 29 15" src="https://github.com/user-attachments/assets/da17496b-2731-4e18-8c35-482236cc63c0" />

  
    
## 07. Try out more

If you like, try out to compare the full genomes of the **[fourspine stickleback](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/048/569/185/GCA_048569185.1_Unibe_ApeQuad_male_1.0/GCA_048569185.1_Unibe_ApeQuad_male_1.0_genomic.fna.gz)** to that of the **[threespine stickleback](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/920/845/GCF_016920845.1_GAculeatus_UGA_version5/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna.gz)**
in Dgenies, zoom into the X and Y of the fourspine stickleback, what do you see?

If you still like to explore more, you can compare the Z chromosomes of the duck (GCA_047663525.1) and ostrich (GCA_040807025.1)
or the full genomes that revealed neo sex chromosome formation in parrots, for this compare the monk parakeet (GCA_017639245.1) to chicken (GCA_024206055.2)

