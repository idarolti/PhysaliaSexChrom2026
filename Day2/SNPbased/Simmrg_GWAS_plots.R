# Plot GWAS

# load required packages

library(tidyverse)

# check working directory
getwd()
# change if necessary

# read in GWAS file, remove mitochondria and unplaced scaffolds
gwas1 <- read_tsv("Simmrg_test/Simmrg3.assoc.txt") %>% filter(grepl("LS",rs))

# split position names to separate chromosomes and point
gwas2 <- mutate(gwas1, scaf = str_split_i(gwas1$rs, ":", 1))

# select necessary columns
gwas <- select(gwas2, "scaf", "ps", "n_mis", "n_obs", "allele0", "allele1", "af", "p_lrt") %>%
  na.omit()

# read in reference index, rename columns to match GWAS and keep only chromosomes
scaffold_lengths <- read_tsv("astcal.fasta.fai", col_names = c("scaf","length")) %>% filter(grepl("LS",scaf))

# set a 5% cut off
# find 5% of the number of values
cutoff5 <- round(nrow(gwas)*0.05)
# sort values by p_lrt, and filter for the first 5%
gwas_sort_cut5 <- head(gwas %>% arrange(p_lrt), n=cutoff5) 
# find the threshold p-value
max(gwas_sort_cut5$`p_lrt`)
perc5 <-max(gwas_sort_cut5$`p_lrt`)

# this value is the 5% cutoff

# do the same for a 1% cutoff
cutoff1 <- round(nrow(gwas)*0.01)
gwas_sort_cut1 <- head(gwas %>% arrange(p_lrt), n=cutoff1) 
max(gwas_sort_cut1$`p_lrt`)
perc1 <-max(gwas_sort_cut1$`p_lrt`)

# the lowest value is the 1% cutoff

# how many values are above the 1% threshold in each chromosome?
table(gwas_sort_cut1$scaf)

# as above but another method
quantile(gwas$`p_lrt`, c(0.95, 0.99), na.rm = T)

# mean p-value
mean(gwas$p_lrt)

# Convert the p-values to the -log10 scale.
gwas4 <- transform(gwas,p_lrt = -log10(p_lrt))

# merge the GWAS values to the reference index. 
gwas_join <- left_join(scaffold_lengths,gwas4) %>% filter(p_lrt>=0)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(gwas_join$scaf))
gwas_join$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(gwas_join$scaf)){
  nbp[i] <- max(gwas_join[gwas_join$scaf == i,]$ps)
  gwas_join[gwas_join$scaf == i,"BPcum"] <- gwas_join[gwas_join$scaf == i,"ps"] + s
  s <- s + nbp[i]
}

axis.set <- gwas_join %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

gwas_join %>% ggplot(aes(x=BPcum,y=p_lrt,color=as.factor(scaf))) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="CHROMOSOME",color="") + 
  guides(color = "none") + 
  geom_hline(yintercept = perc1, linetype="dotted", color = "black", size=0.75) + 
  geom_hline(yintercept = perc5, linetype="dotted", color = "black", size=1.5)  + 
  labs(title="GWAS for P. marginatus (North) with 5% and 1% cutoff")

