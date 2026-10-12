### ------------------------------------------------------
### 0. Download files 
### ------------------------------------------------------

## prerequisites: download the output .fai, .bed and .results file from the SCINKD run to local machine
## depending on the genome used, you might need to get a better naming of the chromosomes


### ------------------------------------------------------
### 1. Prepare files and env
### ------------------------------------------------------

## set directory to the one that has the file downloaded from the server

setwd("~/Desktop/Physalia2026/Day5/scinkd2/scinkd2_AnnSte/")

## install and activate all libraries you will be using
## install.packages("devtools")
library(devtools)
## install from GitHub repository
## devtools::install_github("daewoooo/SVbyEye", branch="master")
library(SVbyEye)
## install.packages("dplyr")
library(dplyr)
## install.packages("ggplot2")
library(ggplot2)
## install.packages("ggrepel")
library(ggrepel)

## inspect chromosome naming and number and potentially fix
## provide number of chromosomes below
## we get this information from the index file ".fai" that was generated for the genome .fasta file with samtools
## that file contains also the length of each entry in the multi fasta file in column 2; for details see here https://www.htslib.org/doc/faidx.html

hap1_idx <- read.delim("RESULTS/rAnnSte1.2_.hap1.fasta.gz.fai", header = F)
hap2_idx <- read.delim("RESULTS/rAnnSte1.2_.hap2.fasta.gz.fai", header = F)

View(hap1_idx)
## in this example, the chromosomes of both haplotypes start with "CM" 
View(hap2_idx)


## subset to major chromosomes
hap1_idx <- hap1_idx[grepl("CM", hap1_idx$V1), ]
## subset to major chromosomes
hap2_idx <- hap2_idx[grepl("CM", hap2_idx$V1), ]

## now we want to add short chromosome names
##In the index file, chromosomes are in the same order as in the fasta and that should correspond to their naming as displayed in the NCBI record
chrnames<-c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9","Chr10")
hap1_idx$chrnames<-chrnames
hap2_idx$chrnames<-chrnames

NCHROM<-10
prefix<-"AnnSte"

### ------------------------------------------------------
### 2. Read kmer files
### ------------------------------------------------------
hap1_kmer <- read.delim("RESULTS/rAnnSte1.2_.hap1-minus-hap2.results", header = F)
hap2_kmer <- read.delim("RESULTS/rAnnSte1.2_.hap2-minus-hap1.results", header = F)

# only get the kmers for the major chromosomes
hap1_tmp <- subset(hap1_idx, V1 %in% hap1_kmer$V1)
hap2_tmp <- subset(hap2_idx, V1 %in% hap2_kmer$V1)

hap1_merge <- merge(hap1_tmp, hap1_kmer, by = "V1", sort = F)
hap2_merge <- merge(hap2_tmp, hap2_kmer, by = "V1", sort = F)

hap1_merge$density <- with(hap1_merge, V2.y / V2.x)
hap2_merge$density <- with(hap2_merge, V2.y / V2.x)


### ------------------------------------------------------
### 3. Subsample first NCHROM entries to plot only the major chromosomes
### ------------------------------------------------------
hap1_final <- hap1_merge[1:NCHROM, c(1,2,6,7,8)]
hap2_final <- hap2_merge[1:NCHROM, c(1,2,6,7,8)]


### ------------------------------------------------------
### 4. Combine & label
### ------------------------------------------------------
hap1_final$Dataset <- "hap1"
hap2_final$Dataset <- "hap2"
combined_data <- rbind(hap1_final, hap2_final)


### ------------------------------------------------------
### 5. Boxplot of kmer densities per haplotype 
### ------------------------------------------------------
## pdf(boxplot_pdf, width = 7, height = 5)
## or
## png(boxplot_png, width = 1600, height = 1200, res = 200)


boxplot(hap1_final$density, hap2_final$density,
        pch = 19,
        names = c("hap1", "hap2"),
        main = paste("Density Comparison (", prefix, ")", sep = ""),
        ylab = "Density")
## dev.off()

## optional: generate Boxplot with outlier labeling

## define outlier function
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
## identifiy outliers
hap1_final$outlier<-ifelse(is_outlier(hap1_final$density), hap1_final$chrnames, NA)
hap2_final$outlier<-ifelse(is_outlier(hap2_final$density), hap2_final$chrnames, NA)

combined_data <- rbind(hap1_final, hap2_final)

## plot
ggplot(combined_data, aes(x = Dataset, y = density)) +
  geom_boxplot() +
  geom_text_repel(aes(label = outlier), na.rm = TRUE) +
  labs(title = (paste("Hap-mer Density ",prefix)),
       x = "Haplotypes",
       y = "Density") +
  theme_minimal()


### ------------------------------------------------------
### 6. Dotplot of kmer density comparing the haplotypes
### ------------------------------------------------------
dotplot <- ggplot(combined_data, aes(x = V2.x, y = V2.y, color = Dataset)) +
  geom_point() +
  geom_smooth(method = "lm", level = 0.95) +
  labs(title = paste("Sequence Length vs # Hap-mers — ", prefix, sep = ""),
       x = "Sequence Length",
       y = "# Hap-mers") +
  theme_minimal()

## pdf(dotplot_pdf, width = 7, height = 5)
## or
## png(dotplot_png, width = 1600, height = 1200, res = 200)
print(dotplot)
## dev.off()


## adjusted plot with point labels
dotplot <- ggplot(combined_data, aes(x = V2.x, y = V2.y, color = Dataset)) +
  geom_point() +
  geom_smooth(method = "lm", level = 0.95) +
  labs(title = paste("Sequence Length vs # Hap-mers — ", prefix, sep = ""),
       x = "Sequence Length",
       y = "# Hap-mers") +
  geom_text(
    label=(combined_data$chrnames), 
    nudge_x = 0.9, nudge_y = 0.9, vjust = -1,    # moves text up
    hjust = 0.5,
    check_overlap = F
  )+
  theme_minimal()

print(dotplot)


### ------------------------------------------------------
### 7. Mirrored kmer density along a specific chromosome
### ------------------------------------------------------

## run below to plot kmers along a specific chromosome
## attention, the input files are quite large and hence loading will take a moment
## you can also get the reduced bed files for the chromosome of interest from the share folder, if you load those, skip the subset command, this would be
## rAnnSte1.2_.Chr10.hap1-minus-hap2.bed
## rAnnSte1.2_.Chr10.hap2-minus-hap1.bed

hap1_bed <- read.delim("RESULTS/rAnnSte1.2_.hap1-minus-hap2.bed", header = F)
hap2_bed <- read.delim("RESULTS/rAnnSte1.2_.hap2-minus-hap1.bed", header = F)
colnames(hap1_bed)<-c("scaffold","start","end")
colnames(hap2_bed)<-c("scaffold","start","end")

## subset provided bed file to chromosome of interest (LGi)
## we saw that Chr10 is an outlier which is CM119782.2 in hap 2 and CM143401.1 in hap 1
LGi<-"Chr10"
hap1_LGi <- subset(hap1_bed, scaffold == "CM143401.1")
hap2_LGi <- subset(hap2_bed, scaffold == "CM119782.2")

## mirror plot
ggplot(hap1_LGi, aes(x=start)) +
  geom_histogram(data = hap1_LGi,  aes(x = start, y =  after_stat(density)), color = "#5F00F6", fill = "black", alpha = 0.8, binwidth=1e6) +
  geom_histogram(data = hap2_LGi, aes(x = start, y = -after_stat(density)), color = "#62CD32", fill = "grey", alpha = 0.8, binwidth=1e6) +
  ggtitle(paste("Kmer Density: ", LGi,prefix, sep=" ")) +
  xlab("Position on chromosome") + 
  ylab("Kmer density") + 
  theme_bw()

### ------------------------------------------------------
### 8. Synteny plot with SVBeye for chromosome of interest
### ------------------------------------------------------

## visualize haplotype alignment of a linkage group of interest

## prerequisite
## align the two haplotypes using minimap2 -> this is the PAF file
## we provide those for you on the server
## minimap2 -x asm5 -t8 -c --eqx --secondary=no Anniella_stebbinsi_HiFi_2024.asm.hic.hap1.fasta.gz Anniella_stebbinsi_HiFi_2024.asm.hic.hap2.fasta.gz > Anniella_stebbinsi_HiFi_2024.asm.hic.paf

## outside of R subset the paf file to alignemnts of your chromosome of interest using awk in your terminal

## awk ' $1 == "CM119782.2" ' rAnnSte1.2.asm.hic.paf  > tmp
## awk ' $6 == "CM143401.1" ' tmp > rAnnSte1.2.asm.hic.iLG.paf
## watch out, the awk command depends on the way you had done the alignment
## to see which chromosomes are in column one you can do cut -f 1 rAnnSte1.2.hic.paf | uniq

## read in PAF
mypaf <- readPaf("../PAFFILES/rAnnSte1.2.asm.hic.iLG.paf")


## plot alignments
plotGenome(paf.table = mypaf)

plotGenome(paf.table = mypaf, color.by = 'identity', perc.identity.breaks = c(85, 90, 95))

plotMiro(paf.table = mypaf, binsize = 100000, genomic.scale = 'Mbp', color.by = 'identity', perc.identity.breaks = c(85, 90, 95))

plotMiro(paf.table = mypaf, color.by = "direction")


### ------------------------------------------------------







