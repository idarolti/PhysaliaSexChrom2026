# SNP density

# load required packages

library(tidyverse)
library(ggpubr)

# check working directory
getwd()
# Here you must be inside the directory with all the .snpden files


# load in the individual SNP densities calculated for each individual
# get the mean male and mean female SNP Density in that window, run permutation tests to see which windows are significantly different.

# load in the male data, by searching for all files that start with "MALE"
myFiles <- list.files(pattern="^MALE*")

# this code reads each file and adds it to a table for all male samples.
# start with a backbone of the first file
# read in the file as a tab-delimited file, including the column names
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
# split the name of the file into chumnks separated by "_"
firstsplit <- strsplit(myFiles[1], "_")
# split the third chunk of this into more chunks, separated by "."
secondsplit <- strsplit(firstsplit[[1]][2], "\\.")
# make a column named with the first chunk of the first split, and the first chunk of the second split
column_ID <- paste(firstsplit[[1]][1],secondsplit[[1]][1],sep="_")
# combine it all together with the column name
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

#loop the other files

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  secondsplit <- strsplit(firstsplit[[1]][2], "\\.")
  column_ID <- paste(firstsplit[[1]][1],secondsplit[[1]][1],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}

# make two colmns from the chromosome name and position
SNPdensity.males <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

# do the same for the females

# load the female data:

myFiles <- list.files(pattern="^FEMALE*")

# build a backbone
backbone <- read_delim(file=myFiles[1],delim = "\t",col_names = T) 
firstsplit <- strsplit(myFiles[1], "_")
secondsplit <- strsplit(firstsplit[[1]][2], "\\.")
column_ID <- paste(firstsplit[[1]][1],secondsplit[[1]][1],sep="_")
backbone.upgrade <- backbone %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)

#loop the other files

for(i in 2:length(myFiles)){
  file <- read_delim(myFiles[i],delim = "\t",col_names = T) 
  firstsplit <- strsplit(myFiles[i], "_")
  secondsplit <- strsplit(firstsplit[[1]][2], "\\.")
  column_ID <- paste(firstsplit[[1]][1],secondsplit[[1]][1],sep="_")
  file.upgrade <- file %>% unite(LOCATION,c("CHROM","BIN_START"),sep=":") %>% dplyr::rename(!!column_ID := "VARIANTS/KB") %>% select(-SNP_COUNT)
  backbone.upgrade <- full_join(backbone.upgrade,file.upgrade,by="LOCATION")
}


SNPdensity.females <- backbone.upgrade %>% separate(LOCATION, c("scaf","base"),sep=":")

#join them
SNPdensity <- full_join(SNPdensity.males,SNPdensity.females)

#make new tibble, change type for base, and replace NAs
SNPdensity.rows <- SNPdensity
SNPdensity.rows$base <- as.numeric(as.character(SNPdensity.rows$base))
SNPdensity.rows <- SNPdensity.rows %>% replace_na(list(Males = 0, Females = 0)) %>% replace(is.na(.), 0)

# use a subsetter to get the male and female average densities 
male_n <- ncol(SNPdensity.rows[ , grepl( "^MALE_" , names( SNPdensity.rows ) ) ])
males_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "^MALE_" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

female_n <- ncol(SNPdensity.rows[ , grepl( "FEMALE_" , names( SNPdensity.rows ) ) ])
females_true <- bind_cols(SNPdensity.rows %>% select(1:2),SNPdensity.rows[ , grepl( "FEMALE_" , names( SNPdensity.rows ) ) ] %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

# merge means together per window, and find the difference
true_SNPdensity <- full_join(males_true,females_true) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)

# summarise
table(true_SNPdensity$mean_MvF_dif > 0)
table(true_SNPdensity$mean_MvF_dif < 0)
table(true_SNPdensity$mean_MvF_dif == 0)

# so, now we would like to permute!
# this means shuffling - we can shuffle the data columns and then use the code above.

# need an initial run of the permutation:

perm <- sample(3:length(SNPdensity.rows))

SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)

males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)

females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)

true_SNPdensity_position_perm <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females)

perm_backbone <- true_SNPdensity_position_perm %>% select(scaf,base,mean_MvF_dif) %>% rename(p1=mean_MvF_dif)

# then we want to run the rest of the loops (since seed wasn't set, your run of this will likely produce a slightly different final file than others but it should not change the overall results)

#this loop takes maybe 5 minutes to run for 1000 permutations.
for(i in 2:1000){
  perm <- sample(3:length(SNPdensity.rows))
  SNPdensity.rows.perm <- SNPdensity.rows %>% select(1:2,perm)
  males_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select(3:(male_n+2)) %>% mutate(mean_Males = rowMeans(.))) %>% select(scaf,base,mean_Males)
  females_perm <- bind_cols(SNPdensity.rows.perm %>% select(1:2),SNPdensity.rows.perm %>% select((male_n+3):(male_n+2+female_n)) %>% mutate(mean_Females = rowMeans(.))) %>% select(scaf,base,mean_Females)
  column_ID <- paste0("p",i)
  perm_upgrade <- full_join(males_perm,females_perm) %>% mutate(mean_MvF_dif = mean_Males - mean_Females) %>% select(scaf,base,mean_MvF_dif) %>% dplyr::rename(!!column_ID := "mean_MvF_dif")
  perm_backbone <- full_join(perm_backbone,perm_upgrade)
}

# join together the permutation p-values with the SNP density values
perm_with_true <- full_join(true_SNPdensity %>% select(scaf,base,mean_MvF_dif),perm_backbone)

# here we are using a conditional mutate where if the MvF difference is > 0 we look for outliers above, if the MvF difference is < 0  we look for outliers below, and if the MvF difference is ==0, we say p value is 1.
perm_with_true_p <- perm_with_true %>% 
  mutate(Pvalue = case_when(mean_MvF_dif > 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] > perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                            mean_MvF_dif < 0 ~ (rowSums(perm_with_true[,grep("p", names(perm_with_true))] < perm_with_true$mean_MvF_dif)+1)/(length(perm_with_true[,grep("p", names(perm_with_true))])+1),
                            mean_MvF_dif == 0 ~ 1)) 

write_tsv(perm_with_true_p,"SNPdensity_perm_with_true_p_Simmrg.txt")

#for step 3, we only need a few columns out of this file and a subset of the genome, so we'll parse it
write_tsv(perm_with_true_p %>% select(scaf,base,mean_MvF_dif,Pvalue), "SNPdensity_SexFindR_Simmrg.txt")
# adding window size for consistency between analyses

# want a look at proportion of scaffold - find the outlier
scaffold_length <- read_tsv("../astcal.fasta.fai",col_names = F) %>% rename(scaf=X1,length=X2)

# filter for significant (p < 0.001) p-values
perm_with_true_p001_snp <- perm_with_true_p %>% filter(Pvalue <= 0.001) 
# count the number of significant windows per chromosome
count_snp_p001_scaffolds <- perm_with_true_p001_snp %>% select(scaf) %>% count(scaf)
# make this value proportional to the size of the chromosome
proportion_count_snp_p001_scaffolds <- left_join(count_snp_p001_scaffolds,scaffold_length) %>% mutate(proportion = (n*10000)/length)

# bar plot
proportion_count_snp_p001_scaffolds %>% 
  ggplot(aes(x=scaf,y=proportion)) + 
  geom_col() + 
  theme(axis.text.x = element_text(angle = 90))

##### SCATTER PLOT

# SNP density
# read in the file you just created, and filter sites for p <= 0.001.
getwd()
SNP <- read_tsv("SNPdensity_SexFindR_Simmrg.txt") %>% 
  filter(Pvalue <= 0.001)

# merge with the reference index
SNPchr <- left_join(scaffold_length, SNP) %>%
  filter(grepl("LS", scaf))

(
  site <-
    SNPchr %>%
    ggplot(aes(x = base / 1000000,
               y = mean_MvF_dif)) +
    geom_point(size = 1) +
    facet_wrap(~scaf, scale = "free") +
    labs(x = "Mb", color = "") +
    labs(y = "10kb SNP density \n(male mean - female mean)", color = "") +
    labs(
      title = bquote("Density of sex-associated SNPs for P. marginatus (North)"),
      subtitle = "Per window",
      color = ""
    ) +
    scale_color_manual(values = cols,
                       guide = "none") +
    scale_y_continuous(
      limits = c(
        min(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25,
        max(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25
      ),
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(size = 10, vjust = .5),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(
        size = 15,
        angle = 90,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      title = element_text(size = 15)
    )
)


# look at interesting chromosomes by filtering the input file by scaffold
(
  site <-
    SNPchr %>%
    filter(scaf == "LS420023") %>%
    ggplot(aes(x = base / 1000000,
               y = mean_MvF_dif)) +
    geom_point(size = 1) +
    facet_wrap(~scaf, scale = "free") +
    labs(x = "Mb", color = "") +
    labs(y = "10kb SNP density \n(male mean - female mean)", color = "") +
    labs(
      title = bquote("Density of sex-associated SNPs for P. marginatus (North)"),
      subtitle = "Per window",
      color = ""
    ) +
    scale_color_manual(values = cols,
                       guide = "none") +
    scale_y_continuous(
      limits = c(
        min(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25,
        max(SNPchr$mean_MvF_dif, na.rm = TRUE) * 1.25
      ),
      expand = c(0, 0)
    ) +
    theme_bw() +
    theme(legend.position = "none") +
    theme(
      axis.text.x = element_text(size = 10, vjust = .5),
      axis.title.x = element_text(size = 15),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(
        size = 15,
        angle = 90,
        hjust = .5,
        vjust = .5,
        face = "plain"
      ),
      title = element_text(size = 15)
    )
)
