# Plot Fst

# load required packages

library(tidyverse)

# check working directory
getwd()
# change if necessary

# load you Fst ouput
Fst <- read_tsv("Simmrg_test/Simmrg_astcal.windowed.weir.fst") %>% 
  replace_na(list(WEIGHTED_FST=0)) %>% # make any NA values 0
  rename(scaf = CHROM, start = BIN_START, end = BIN_END) %>% # rename files to match to reference index
  mutate(start=as.numeric(start)) %>% # make positions numeric
  filter(WEIGHTED_FST >=0) # discard negative values

# load reference index, rename columns to match fst input, and discard mitochrondria and unplaced scaffolds
scaffold_lengths <- read_tsv("astcal.fasta.fai", col_names = c("scaf","length")) %>% filter(grepl("LS",scaf))

# set a 5% cut off
# find 5% of the number of values
cutoff1 <- round(nrow(Fst)*0.05)
# sort values by Fst, and filter for the first 5%
Fst_sort_cut1 <- head(Fst %>% arrange(-WEIGHTED_FST), n=cutoff1) 
# find the Fst value of this 5% threshold
min(Fst_sort_cut1$`WEIGHTED_FST`)
perc5 <- min(Fst_sort_cut1$`WEIGHTED_FST`)

# the lowest value is the 5% cutoff

# set a 1% cutoff in the same way
cutoff2 <- round(nrow(Fst)*0.01)
Fst_sort_cut2 <- head(Fst %>% arrange(-WEIGHTED_FST), n=cutoff2) 
min(Fst_sort_cut2$`WEIGHTED_FST`)
perc1 <- min(Fst_sort_cut2$`WEIGHTED_FST`)

# the lowest value is the 1% cutoff

# how many values are above the 1% threshold in each chromosome?
table(Fst_sort_cut2$scaf)

# as above but another method
quantile(Fst$`WEIGHTED_FST`, c(0.95, 0.99), na.rm = T)

# mean fst
mean(Fst$WEIGHTED_FST)

# merge the fst values to the reference index. 
# left_join, because we want to keep all rows in the scaffold_lengths, even if they don't appear in Fst, but discard rows in Fst that don't appear in scaffold_length. This way, we remove unplaced scaffolds and mitochondrial positions.
Fst_join <- left_join(scaffold_lengths,Fst)

#for x axis, we want cumulative bases for each position in the genome for a continuous axis
nCHR <- length(unique(Fst_join$scaf))
Fst_join$BPcum <- 0
s <- 0
nbp <- c()
for (i in unique(Fst_join$scaf)){
  nbp[i] <- max(Fst_join[Fst_join$scaf == i,]$start)
  Fst_join[Fst_join$scaf == i,"BPcum"] <- Fst_join[Fst_join$scaf == i,"start"] + s
  s <- s + nbp[i]
}

# give positions for the axes, based on the size of each chromosome
axis.set <- Fst_join %>% 
  group_by(scaf) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)

# plot in ggplot
Fst_join %>% ggplot(aes(x=BPcum,y=WEIGHTED_FST,color=as.factor(scaf))) +
  geom_point(alpha = 0.75) +
  scale_x_continuous(label = axis.set$scaf, breaks = axis.set$center) + 
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x="CHROMOSOME",color="") + 
  guides(color = "none") + 
  geom_hline(yintercept = perc1, linetype="dotted", color = "black", size=0.75) + 
  geom_hline(yintercept = perc5, linetype="dotted", color = "black", size=1.5)  + 
  labs(title="Fst for P. marginatus (North) with 5% and 1% cutoff")

