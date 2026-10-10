# Day 3 Practical - 02. Finding gametolog sequences

This part of the practical will cover the steps for identifying sex-linked sequences with **[SEX-DETector](https://pmc.ncbi.nlm.nih.gov/articles/PMC5010906/)**

## 00. Prepare work folder for day 3

```
ssh -i ~/YOURLOCALFOLDER/scverse1.pem ubuntu@35.89.239.29
mkdir day3
cd day3
conda activate /opt/conda-envs/day3
```

## 02. Genotyping

First, obtain alignment bam files and transcriptome assembly.
```
mkdir sexdetector
cd sexdetector
cp -r ~/Share/day3/sexdetector/transcriptome_assembly/ ./
cp -r ~/Share/day3/sexdetector/bam_files/ ./
cp -r ~/Share/day3/sexdetector/scripts/ ./
```

Genotyping will be done using **[reads2snp](https://kimura.univ-montp2.fr/popphyl/resources/tools/)**. This is preferred when using SEX-DETector with RNA-seq data because it allows for allelic expression biases (important for sex chromosome studies because Y or W alleles may be less expressed).

Make a list of bam files to run reads2snp.

```
cd bam_files
ls -d "$PWD"/* > bam_list.txt
```

Run reads2snp:
- -aeb (allows alleles to have different expression levels, which is important for sex chromosome anaylses as the Y copy can be less expressed than the X copy)
- -min: 3 (minimum number of reads to call a genotype)
- -par: 0 (do not remove SNPs that appear to come from paralogous sequences, avoid overfiltering as X/Y SNPs can look like paralogous SNPs)
- -bqt: 20 (minimum base quality)
- -rqt: 10 (minimum read mapping quality

```
cd ../scripts
mkdir ../reads2snp

./reads2snp_2.0.64.bin <PARAMETERS> -bamlist </PATH/BAM_LIST> -bamref </PATH/ASSEMBLY>.fasta -out ../reads2snp/reads2snp_output
```

<details>
<summary>🔑 View Full Code</summary>

```
./reads2snp_2.0.64.bin -aeb -min 3 -par 0 -bqt 20 -rqt 10 -bamlist ../bam_files/bam_list.txt -bamref ../transcriptome_assembly/trinity.fasta -out ../reads2snp/reads2snp_output
```
</details>

Have a look at the two main reads2snp outputs: .alr and .gen

```
cd ../reads2snp
head reads2snp_output.alr
head reads2snp_output.gen
```

The output .alr file is a tab delimited file that shows each position of each gene on a different line with the major allele (most common allele in number of reads), then M if the position is monomorphic or P if the position shows different alleles, then for each individual the total read number if the position is monomorphic or the total read number followed by the number of reads for each base [A/C/G/T] if the position is polymorphic:

```
>contig_name
	maj	M/P	individual1	individual2
	T	M	8	1
	C	P	17[0/0/0/17]	14[0/6/0/8]
```

The output .gen file is a tab delimited file showing each position of each gene on a separate line with the position number starting from 1, then for each individual the inferred genotype followed by a pipe and the posterior probability of the inferred genotype:

```
>contig_name
  position	individual1	individual2
  1	TT|1	NN|0
  2	TT|1	TT|0.98
```

## 03. SNP segregation analysis (subset)

Identify sex-linked sequences with SEX-DETector. The software can be found **[here](https://gitlab.in2p3.fr/sex-det-family)**. For this practical we will use the original version of the software. Newer versions can also identify the type of sex chromosome system (XY or ZW) if you don't know it a priori.

```
mkdir ../sexdetector_output
cd ../scripts
```

Generate gen_summary file - allows SEX-DETector to run faster. The following options should be use for the command below, -hom string (names of homogametic progeny individuals separated by commas), -het string (the names of the heterogametic progeny individuals separated by commas), -hom_par string (homogametic parent name), -het_par string (heterogametic parent name). The gen_summary file is similar to the gen file except that it shows only one occurence of each possible SNP in the dataset and shows the number of times it happens in the first column instead of the position number of the SNP.

```
./SEX-DETector_prepare_file.pl </PATH/FILE>.gen ../reads2snp/reads2snp_output.gen_summary -hom <NAME_HOMOGAMETIC_INDIVIDUALS> -het <NAME_HETEROGAMETIC_INDIVIDUALS> -hom_par <NAME_HOMOGAMETIC_PARENT> -het_par <NAME_HETEROGAMETIC_PARENT>

head ../reads2snp/reads2snp_output.gen_summary
```

<details>
<summary>🔑 View Full Code</summary>

```
./SEX-DETector_prepare_file.pl ../reads2snp/reads2snp_output.gen ../reads2snp/reads2snp_output.gen_summary -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father
```
</details>


Run SEX-DETector

```
./SEX-DETector.pl -help

./SEX-DETector.pl -alr </PATH/FILE>.alr -alr_gen </PATH/FILE>.gen -alr_gen_sum </PATH/FILE>.gen_summary -system xy -hom <NAME_HOMOGAMETIC_INDIVIDUALS> -het <NAME_HETEROGAMETIC_INDIVIDUALS> -hom_par <NAME_HOMOGAMETIC_PARENT> -het_par <NAME_HETEROGAMETIC_PARENT> -seq -detail -detail-sex-linked -out ../sexdetector_output/Poecilia_reticulata

```

<details>
<summary>🔑 View Full Code</summary>

```
./SEX-DETector.pl -alr ../reads2snp/reads2snp_output.alr -alr_gen ../reads2snp/reads2snp_output.gen -alr_gen_sum ../reads2snp/reads2snp_output.gen_summary -system xy -hom Female_Offspring1,Female_Offspring2,Female_Offspring3,Female_Offspring4,Female_Offspring5 -het Male_Offspring1,Male_Offspring2,Male_Offspring3,Male_Offspring4,Male_Offspring5 -hom_par Female_Mother -het_par Male_Father -seq -detail -detail-sex-linked -out ../sexdetector_output/Poecilia_reticulata
```
</details>

The main outputs from SEX-DETector are:
- SNPs_detail.txt: information on each SNP of each gene (position in gene, likelihood autosomal/sex-linked, genotypes)
- sex-linked_detail.txt: information on SNP inferred as sex-linked (position, likelihood, expression of each allele, error rates)
- sex-linked_sequences.fasta: X and Y sequences fasta file
- assignment.txt: outputs the probability for each gene to be autosomal, sex-linked (with X/Y alleles), sex-linked (XO hemizygous)
  
Look at the outputs produced by SEX-DETector. How many genes are inferred to be sex-linked?


## 04. SNP segregation analysis (full dataset)

Next, get the SEX-DETector output for the full dataset.

```
cd ~/day3/sexdetector
cp -r ~/Share/day3/sexdetector/sexdetector_output_full_RetFam1 ./
cd sexdetector_output_full_RetFam1
```

Find how many genes are inferred as autosomal versus sex-linked.

```
cat RetFam1_assignment.txt | grep -w "sex-linked" -c
cat RetFam1_assignment.txt | grep -w "autosomal" -c
```

Obtain a single sequence for each gene (to then blast onto the assembly).

```
head RetFam1_sex-linked_sequences.fasta

# Reads file line by line
awk '
	# For lines that start with >
	/^>/ {
		# Save the header line in a variable
		gene = $0
		# Remove the leading > sign
		sub(/^>/, "", gene)
		# Remove the suffix after the underscore
		sub(/_[^_]+$/, "", gene)
		# Use an array to check if this Gene ID was seen before
		if (seen[gene]++) {
			# If the Gene ID has been seen already then skip
			skip = 1
		# If the Gene ID has not been seen
		} else {
			# Print the header
			print $0
			# Skip becomes 0
			skip = 0
		}
	}
	# For lines that do not start with >, the sequence lines
	!/^>/ {
		# If skip is 0, print sequence line
		if (!skip) print $0
		}
' RetFam1_sex-linked_sequences.fasta > RetFam1_sex-linked_sequences_unique.fasta
			
grep ">" RetFam1_sex-linked_sequences_unique.fasta -c
```

Use Blast to see where on the genome the identifyied sex-linked genes align. A Blast database (with makeblastdb) has been already created in the Shared folder.

```
blastn -db ~/Share/day3/sexdetector/genome_assembly/Poecilia_reticulata -query RetFam1_sex-linked_sequences_unique.fasta -out RetFam1_blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

head RetFam1_blastout
```

Identify top blast hits for each sequence. This script takes a blast output file (format: outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq") and identifies the top blast hit for each query. Top blast hit = minimum 30 pidentity, greatest blast score and greatest pidentity. If a query has two hits with identical blast score and pidentity, one is chosen randomly as the tophit. 

```
# Keep only hits with percent identity > 30 (column 3)
awk '$3 > 30' RetFam1_blastout \
  | \
  # Sort by Gene ID (column 1), then Bitscore descending (column 12), then % Identity descending (column 3)
  sort -k1,1 -k12,12nr -k3,3nr \
  | \
  # Pick the top hit per gene (skip if hit #1 and hit #2 tie)
  awk -v OFS=',' '
    function print_top() {
        if (gene != "" && !is_tied) {
            print gene, subject, bitscore, pident, sstart, send
        }
    }
    # When we hit a new Gene ID
    $1 != gene {
        print_top()
        gene=$1; subject=$2; pident=$3; sstart=$9; send=$10; bitscore=$12
        is_tied = 0
        next
    }
    # If the second hit has the exact same bitscore and identity, mark as tied
    $12 == bitscore && $3 == pident {
        is_tied = 1
    }
    END { print_top() }
  ' > RetFam1_blastout_tophits

head RetFam1_blastout_tophits
```

Count the number of times a chromosome is assigned a sex-linked gene.

```
# Read file line by line, specifying the files is comma-separated
awk -F',' '
{
    # Count occurrences of each chromosome (column 2)
    genes_per_chrom[$2]++
}
END {
    for (chrom in genes_per_chrom) {
        print chrom, genes_per_chrom[chrom]
    }
}
' RetFam1_blastout_tophits
```

Extract the inferred sex-linked genes that align to the sex chromosome (CM002717.1).

```
awk -F',' '$2 == "CM002717.1"' RetFam1_blastout_tophits > RetFam1_blastout_tophits_sexchromo
```

Transfer this file to your local machine, and plot the distribution of sex-linked genes across the sex chromosome using R.

```
library(ggplot2)

# Load data and assign column names
sexlinked <- read.csv("RetFam1_blastout_tophits_sexchromo", header = FALSE)
names(sexlinked) <- c("Gene", "Chromosome", "Bitscore", "PIdentity", "Start", "End")

# Density plot
ggplot(sexlinked, aes(x = Start / 1e6)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  # Add tick marks for the position of sex-linked genes
  geom_rug(color = "firebrick", length = unit(0.06, "npc"), linewidth = 0.8) +
  xlim(0, 26) +
  labs(
    title = "Sex-Linked Gene Distribution",
    subtitle = paste0("Chromosome: ", sexlinked$Chromosome[1], " (N = ", nrow(sexlinked), " genes)"),
    x = "Chromosome Position (Mb)",
    y = "Density"
  ) +
  theme_classic()
```


**Task: Run the last part of the analysis (04.SNP segregation analysis full) using the SEX-DETector output from another family.**

For comparison, assess the quality of your cleaned data:
1. Copy to your directory the folder ~/Share/day3/sexdetector/sexdetector_output_full_RetFam2
2. Find how many genes are inferred as autosomal versus sex-linked
3. Obtain a single (unique) sequence for each gene
4. Blast unique sequences to the assembly and identify top blast hits
5. Extract "true" sex-linked genes (as those aligning to the sex chromosome CM002717.1)
6. Transfer the output to your desktop and use R to plot the distribution of sex-linked genes across the chromosome

**What differences do you find in the number and distribution of sex-linked genes between Family 1 and Family 2?**

<details>
<summary>🔑 View Solution</summary>

```
cd ~/day3/sexdetector/
cp -r ~/Share/day3/sexdetector/sexdetector_output_full_RetFam2 ./
cd sexdetector_output_full_RetFam2
cat RetFam2_assignment.txt | grep -w "sex-linked" -c
cat RetFam2_assignment.txt | grep -w "autosomal" -c

# Reads file line by line
awk '
	# For lines that start with >
	/^>/ {
		# Save the header line in a variable
		gene = $0
		# Remove the leading > sign
		sub(/^>/, "", gene)
		# Remove the suffix after the underscore
		sub(/_[^_]+$/, "", gene)
		# Use an array to check if this Gene ID was seen before
		if (seen[gene]++) {
			# If the Gene ID has been seen already then skip
			skip = 1
		# If the Gene ID has not been seen
		} else {
			# Print the header
			print $0
			# Skip becomes 0
			skip = 0
		}
	}
	# For lines that do not start with >, the sequence lines
	!/^>/ {
		# If skip is 0, print sequence line
		if (!skip) print $0
		}
' RetFam2_sex-linked_sequences.fasta > RetFam2_sex-linked_sequences_unique.fasta

grep ">" RetFam2_sex-linked_sequences_unique.fasta -c

blastn -db ~/Share/day3/sexdetector/genome_assembly/Poecilia_reticulata -query RetFam2_sex-linked_sequences_unique.fasta -out RetFam2_blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq"

# Keep only hits with percent identity > 30 (column 3)
awk '$3 > 30' RetFam2_blastout \
  | \
  # Sort by Gene ID (column 1), then Bitscore descending (column 12), then % Identity descending (column 3)
  sort -k1,1 -k12,12nr -k3,3nr \
  | \
  # Pick the top hit per gene (skip if hit #1 and hit #2 tie)
  awk -v OFS=',' '
    function print_top() {
        if (gene != "" && !is_tied) {
            print gene, subject, bitscore, pident, sstart, send
        }
    }
    # When we hit a new Gene ID
    $1 != gene {
        print_top()
        gene=$1; subject=$2; pident=$3; sstart=$9; send=$10; bitscore=$12
        is_tied = 0
        next
    }
    # If the second hit has the exact same bitscore and identity, mark as tied
    $12 == bitscore && $3 == pident {
        is_tied = 1
    }
    END { print_top() }
  ' > RetFam2_blastout_tophits

# Read file line by line, specifying the files is comma-separated
awk -F',' '
{
    # Count occurrences of each chromosome (column 2)
    genes_per_chrom[$2]++
}
END {
    for (chrom in genes_per_chrom) {
        print chrom, genes_per_chrom[chrom]
    }
}
' RetFam2_blastout_tophits

awk -F',' '$2 == "CM002717.1"' RetFam2_blastout_tophits > RetFam2_blastout_tophits_sexchromo

# transfer to local directory and run in R

library(ggplot2)

# Load data and assign column names
sexlinked <- read.csv("RetFam2_blastout_tophits_sexchromo", header = FALSE)
names(sexlinked) <- c("Gene", "Chromosome", "Bitscore", "PIdentity", "Start", "End")

# Density plot
ggplot(sexlinked, aes(x = Start / 1e6)) +
  geom_density(fill = "steelblue", alpha = 0.4) +
  # Add tick marks for the position of sex-linked genes
  geom_rug(color = "firebrick", length = unit(0.06, "npc"), linewidth = 0.8) +
  xlim(0, 26) +
  labs(
    title = "Sex-Linked Gene Distribution",
    subtitle = paste0("Chromosome: ", sexlinked$Chromosome[1], " (N = ", nrow(sexlinked), " genes)"),
    x = "Chromosome Position (Mb)",
    y = "Density"
  ) +
  theme_classic()
```
</details>

