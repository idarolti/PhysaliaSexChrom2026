# Day 1 Practical

This practical will cover:

1. Read quality check and trimming
2. Read mapping to genome assembly
3. Variant calling

## 00. Prepare work folder for day 1
```
ssh -i scverse1.pem ubuntu@35.89.239.29
mkdir day1
cd day1
conda activate /opt/conda-envs/day1
```

## 01. Read quality check and trimming

* **[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** - A quality control tool for high throughput sequence data

First, create an output directory for FastQC. Then, run FastQC for all the paired fastq files.

```
mkdir 01.trimming
cd 01.trimming
mkdir fastqc_output_raw_reads
cd fastqc_output_raw_reads

for f in /home/ubuntu/Share/day1/01.quality_trimming/raw_reads/*fastq; do fastqc $f --outdir ./; done
```

* **[MultiQC](https://multiqc.info)** - A tool for merging FastQC output reports of individual samples into a single summary report

This software uses as input the fastqc.zip files produced by FastQC. After running, download the .html output file to your local machine to visualize the results in a web browser.

```
multiqc ./ -o ./
```

After running, download the .html output file to your local machine to visualize the results in a web browser.

```
pwd
scp -i scverse1.pem ubuntu@35.89.239.29:/pwd/*.html ~/Desktop
```

Can find nice examples of different fastqc outputs [here](https://rtsf.natsci.msu.edu/genomics/technical-documents/fastqc-tutorial-and-faq.aspx).

* **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)** - A read trimming tool for Illumina NGS data

The following command will trim reads to remove adapter sequences, regions where the average Phred score in sliding windows of four bases is <15, reads for which the leading/trailing bases have a Phred score <3, and paired-end reads where either read pair is <50 bp. You can find adapter sequences [here](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/SHARE/FrontPages/AdapterSeq.htm). You can use the command on each pair (forward/R1 + reverse/R2) of fastq files.

```
cd ../
mkdir trimmed_reads
    
input_dir="/home/ubuntu/Share/day1/01.quality_trimming/raw_reads/"
adapter_dir="/home/ubuntu/Share/day1/01.quality_trimming/"
output_dir=./trimmed_reads
    
trimmomatic PE \
   $input_dir/sample1_R1.fastq \
   $input_dir/sample1_R2.fastq \
   $output_dir/sample1_output_R1_paired.fastq.gz \
   $output_dir/sample1_output_R1_unpaired.fastq.gz \
   $output_dir/sample1_output_R2_paired.fastq.gz \
   $output_dir/sample1_output_R2_unpaired.fastq.gz \
   ILLUMINACLIP:$adapter_dir/adapters.fa:2:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50
```

**Task: Run FastQC and MultiQC on trimmed reads**

For comparison, assess the quality of your cleaned data:
1. Create a directory named `fastqc_output_trimmed_reads`.
2. Run `fastqc` on all paired files inside `./trimmed_reads/` and direct outputs to your new folder.
3. Run `multiqc` to aggregate the resulting reports into the same folder.

<details>
<summary>💡 Need a hint?</summary>

- Use `mkdir` to create the destination folder.
- In bash, you can iterate over files with:  
  `for f in ./path/to/files/*.fastq.gz; do <command>; done`  
- Both `fastqc` and `multiqc` use `--outdir` to designate the output directory.

</details>

<details>
<summary>🔑 View Solution</summary>

```
mkdir fastqc_output_trimmed_reads
for f in ./trimmed_reads/*_paired.fastq.gz; do fastqc "$f" --outdir ./fastqc_output_trimmed_reads; done
multiqc ./fastqc_output_trimmed_reads --outdir ./fastqc_output_trimmed_reads
```

</details>

Note - for RAD-seq data, trimming also includes removing barcodes and restriction site remnants. Tools like _process_radtags_ from **[STACKS](https://catchenlab.life.illinois.edu/stacks/)** can specifically handle these steps.

## 02. Read mapping

* **[Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)** - A tool for aligning short-read data to a reference genome or genomic sequences

First, build an index for the reference genome. This takes a while to run (SKIP)!

```
bowtie2-build Poecilia_picta.fna Poecilia_picta
```

Then, align each pair of reads to the indexed genome using bowtie2 and convert the output alignment sam file into a sorted bam file.

```
cd ../
mkdir 02.read_alignments
cd 02.read_alignments

bowtie2 -p4 -x /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta \
   -1 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R1_subset.fastq \
   -2 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_female1_R2_subset.fastq \
   | samtools view -b -S - | samtools sort - -o ./Poecilia_picta_female1_subset.bam
```

Filter alignment files by mapping quality.

```
samtools view -b -q 30 Poecilia_picta_female1_subset.bam > Poecilia_picta_female1_subset_mapq.bam

samtools flagstat Poecilia_picta_female1_subset.bam
samtools flagstat Poecilia_picta_female1_subset_mapq.bam
```

Filter alignment files by uniquely mapping reads.

Bowtie2 XS flag is used for this, though it can be ambiguous and inconsistent! 

Earlier BWA version (aln/sampe) XT:A:U flag is more reliable.

HISAT2 also has useful options to filter alignments (--no-discordant --no-mixed --no-unal)

**Task: Align reads and filter alignment file for male sample**

1. Align reads with bowtie2
2. Convert SAM to BAM and sort with samtools
3. Filter alignment for a minimum mapping quality of 30 with samtools
   
- The fastq files for the male sample (Poecilia_picta_male2) are also located in /home/ubuntu/Share/day1/02.read_mapping/reads/

<details>
<summary>💡 Need a hint?</summary>

- bowtie2 takes as input the indexed genome with -x option, the forward and reverse reads with -1 and, respectively, -2 options
- the full path to all input files must be provided, not just the file names
- you can pipe (|) the output alignment file through samtools to convert SAM to BAM and sort the alignments
- use samtools view with -q 30 option to filter alignment file

</details>

<details>
<summary>🔑 View Solution</summary>

```
bowtie2 -p4 -x /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta \
   -1 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_male2_R1_subset.fastq \
   -2 /home/ubuntu/Share/day1/02.read_mapping/reads/Poecilia_picta_male2_R2_subset.fastq \
   | samtools view -b -S - | samtools sort - -o ./Poecilia_picta_male2_subset.bam

samtools view -b -q 30 Poecilia_picta_male2_subset.bam > Poecilia_picta_male2_subset_mapq.bam
```

</details>

## Visualize alignments

Install **[IGV](https://igv.org/doc/desktop/#DownloadPage/)** locally.

Convert .bam file to .bw format, which allows for easier visualization of number of reads mapping in each genomic region. 
(If these steps take too long to run, you can copy the bw files directly).

Bam files for a single chromosome can be found at /home/ubuntu/Share/day1/02.read_mapping/read_alignments/ (Poecilia_picta_female1_chr12.bam, Poecilia_picta_male2_chr12.bam).

**Task: Convert to bw format for visualization**

- **[bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html)**: specify input bam file with -b option, number of threads to run with -p option (can choose 8), and output bw file with -o option
- full path to files is required
- once you have the bw files, transfer them to desktop, together with the genome fasta and index files (found at /home/ubuntu/Share/day1/02.read_mapping/reference_genome/)
- download also the full genome bw files (Poecilia_picta_female1.bw and Poecilia_picta_male2.bw), which can be found at /home/ubuntu/Share/day1/02.read_mapping/read_alignments/

<details>
<summary>🔑 View Solution</summary>

```
bamCoverage -p 8 -b /home/ubuntu/Share/day1/02.read_mapping/read_alignments/Poecilia_picta_female1_chr12.bam -o Poecilia_picta_female1_chr12.bw

scp -i chrsex25.pem ubuntu@44.249.25.243:/home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna ~/Desktop
scp -i chrsex25.pem ubuntu@44.249.25.243:/home/ubuntu/Share/day1/02.read_mapping/read_alignments/Poecilia_picta_female1.bw ~/Desktop
scp -i chrsex25.pem ubuntu@44.249.25.243:/home/ubuntu/Share/day1/02.read_mapping/read_alignments/Poecilia_picta_male2.bw ~/Desktop
```

</details>


Open IGV and load the genome and .bw files.

In IGV, change track height to 200, change color, select type graph "Bar Chart", Window function "Mean", select log scale.

Look at the read mapping rates for the different chromosomes. Can you identify the sex chromosome?



## 03. Variant calling

* **[GATK](https://gatk.broadinstitute.org/hc/en-us)** - A genomic analysis toolkit focused on variant discovery.

Create sequence dictionary for the reference sequence. This takes a while to run (SKIP)!

```
gatk CreateSequenceDictionary -R Poecilia_picta.fna -O Poecilia_picta.dict
```

Call SNPs using [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller). This first part takes a while to run, so we will only use the subset alignment files.

```
mkdir 03.snp_calling

picard AddOrReplaceReadGroups I=./02.read_alignments/Poecilia_picta_female1_subset_mapq.bam O=./02.read_alignments/Poecilia_picta_female1_subset_mapq_RG.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=picta_female1

samtools index ./02.read_alignments/Poecilia_picta_female1_subset_mapq_RG.bam

gatk HaplotypeCaller \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -I ./02.read_alignments/Poecilia_picta_female1_subset_mapq_RG.bam \
   -O ./03.snp_calling/Poecilia_picta_female1_subset.gvcf --emit-ref-confidence GVCF \
   --min-base-quality-score 30 --pcr-indel-model NONE --sample-name picta_female1
```

Perform genotyping of variants using [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/13832766863259-GenotypeGVCFs). The next steps run more quickly, so we can use as input file a gvcf based on an entire chromosome.

```
cp /home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_female1_chr12.gvcf ./03.snp_calling
cp /home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_female1_chr12.gvcf.idx ./03.snp_calling

gatk GenotypeGVCFs \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   --variant ./03.snp_calling/Poecilia_picta_female1_chr12.gvcf \
   -O ./03.snp_calling/Poecilia_picta_female1_chr12.genotyped.gvcf
```

Filter variants using [SelectVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants) and [VariantFiltration](https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration)

```
gatk SelectVariants \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./03.snp_calling/Poecilia_picta_female1_chr12.genotyped.gvcf \
   -O ./03.snp_calling/Poecilia_picta_female1_chr12.selectvar.gvcf --restrict-alleles-to BIALLELIC --select-type-to-include SNP

gatk VariantFiltration \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./03.snp_calling/Poecilia_picta_female1_chr12.selectvar.gvcf \
   -O ./03.snp_calling/Poecilia_picta_female1_chr12.selectvar_filtered.gvcf \
   --filter-expression "QUAL <= 30.0 || DP <= 20" --filter-name "low_qual_or_dp"
```

Obtain the number of homozygous and heterozygous sites

```
vcftools --vcf ./03.snp_calling/Poecilia_picta_female1_chr12.selectvar_filtered.gvcf --extract-FORMAT-info GT --stdout | awk 'NR > 1 {print $3}' | tr '|' '/' | sort | uniq -c
```

**Task: Run genotyping and filtering for the male sample**

1. Copy the male gvcf file (/home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_male2_chr12.gvcf) to your folder
2. Run variant genotyping with gatk GenotypeGVCFs
3. Filter variants with gatk SelectVariants and VariantFiltration
4. Obtain the number of homozygous and heterozygous sites

<details>
<summary>🔑 View Solution</summary>

```
cp /home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_male2_chr12.gvcf ./03.snp_calling
cp /home/ubuntu/Share/day1/03.snp_calling/Poecilia_picta_male2_chr12.gvcf.idx ./03.snp_calling

gatk GenotypeGVCFs \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   --variant ./03.snp_calling/Poecilia_picta_male2_chr12.gvcf \
   -O ./03.snp_calling/Poecilia_picta_male2_chr12.genotyped.gvcf

gatk SelectVariants \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./03.snp_calling/Poecilia_picta_male2_chr12.genotyped.gvcf \
   -O ./03.snp_calling/Poecilia_picta_male2_chr12.selectvar.gvcf --restrict-alleles-to BIALLELIC --select-type-to-include SNP

gatk VariantFiltration \
   -R /home/ubuntu/Share/day1/02.read_mapping/reference_genome/Poecilia_picta.fna \
   -V ./03.snp_calling/Poecilia_picta_male2_chr12.selectvar.gvcf \
   -O ./03.snp_calling/Poecilia_picta_male2_chr12.selectvar_filtered.gvcf \
   --filter-expression "QUAL <= 30.0 || DP <= 20" --filter-name "low_qual_or_dp"

vcftools --vcf ./03.snp_calling/Poecilia_picta_male2_chr12.selectvar_filtered.gvcf --extract-FORMAT-info GT --stdout | awk 'NR > 1 {print $3}' | tr '|' '/' | sort | uniq -c
```

</details>
