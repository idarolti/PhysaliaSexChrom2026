### ------------------------------------------------------
### 0. Download files 
### ------------------------------------------------------

## prerequisites: download the output .fai, .bed and .results file from the SCINKD run to local machine
## depending on the genome used, you might need to get a better naming of the chromosomes


### ------------------------------------------------------
### 1. Prepare files and env
### ------------------------------------------------------

## set directory to the one that has the file downloaded from the server

setwd("~/Desktop/Physalia2026/Day5/scinkd2/scinkd2_SphNot/")

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

hap1_idx <- read.delim("RESULTS/SphNot.hap1.fasta.gz.fai", header = F)
hap2_idx <- read.delim("RESULTS/SphNot.hap2.fasta.gz.fai", header = F)

View(hap1_idx)
## in this example, the names don't give away which of the scaffolds could be chromosomes
## a google search tells that this species has 17 chromosomes so we reduce the file to that number
View(hap2_idx)

## subset to major chromosomes
hap1_idx <- hap1_idx[1:17, ]
## subset to major chromosomes
hap2_idx <- hap2_idx[1:17, ]  

## now we want to add short chromosome names
## in the index file, chromosomes are in the same order as in the fasta and that should correspond to their naming
hap1_idx$chrnames <- paste0("Chr", rownames(hap1_idx))
hap2_idx$chrnames <- paste0("Chr", rownames(hap2_idx))

NCHROM <- 17
prefix <- "SphNot"

### ------------------------------------------------------
### 2. Read kmer files
### ------------------------------------------------------
hap1_kmer <- read.delim("RESULTS/SphNot.hap1-minus-hap2.results", header = F)
hap2_kmer <- read.delim("RESULTS/SphNot.hap2-minus-hap1.results", header = F)

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

## we already see that the two haplotypes are very different in kmer density
## this is suggesting issues with the assembly


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

##the slopes again differ, there is a difference between the two haplotypes

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








