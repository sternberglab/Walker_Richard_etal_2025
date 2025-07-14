
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggtree)
library(biomartr)
library(DECIPHER)

###################################################################################
# Check for reads that show evidence of tail fiber switching in BIDMC93 (sSL3690)

setwd("/Users/tannerwiegand/Documents/Projects/TnpB/dTnpB/Flagellin_candidate_selection/Strains_we_can_get/Mutant_with_prophage_tail_inversion_WGS/HTS_reads")
###################################################################################

# Start clean
rm(list = ls())

# Import FASTQ from plasmidsaurus
# reads <- readDNAStringSet("../Sternberg_S2F_1/reads/raw_reads.fastq", format = "fastq")
#reads <- readDNAStringSet("SRR2127654.fastq", format = "fastq")
#reads <- readDNAStringSet("SRR2135242.fastq", format = "fastq")
reads <- readDNAStringSet("Sternberg_X6Z_raw_reads/Sternberg_X6Z_1_vSL089_BigLinear_53.fastq", format = "fastq")

  
# Add the different fragments that you'll use to search
static <- "GACCGGGGCAACGAATGTCGCTGACGCTCGCACAAACCTC"
tfp1 <- "GGTTTAGGAACATCAGCCATACTTAATGCGCGGTCCAACG"
tfp2 <- "GGTTTGGTAGACAGCAATGGTTACGTGCCTGTGTCACTGG"
tfp3 <- "GGTTTAGGAAGTAGCGCGACACGGGATGCTTACAGCTCGA"
tfp4 <- "GACCTTTATTCACCCGCATCAGCAGTCATGGCCAGCTCTG"


# Try shortening the different fragments that you'll use to search
total.length <- 25 # length upstream and downstream of junction
static <- str_sub(static, -total.length, -1)
tfp1 <- str_sub(tfp1, 1, total.length)
tfp2 <- str_sub(tfp2, 1, total.length)
tfp3 <- str_sub(tfp3, 1, total.length)
tfp4 <- str_sub(tfp4, 1, total.length)


#### Check for TFP1
x <- vmatchPattern(paste0(static,tfp1), reads, max.mismatch = 5)
tfp1.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp1.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp1.df$id])))
tfp1.df <- tfp1.df %>% left_join(y, by = "id")
tfp1.df$match.seq <- NA
tfp1.df$nmismatch <- NA
if(nrow(tfp1.df) > 0){
for(i in 1:nrow(tfp1.df)){ 
  tfp1.df$match.seq[i] <- str_sub(tfp1.df$seq[i], tfp1.df$start[i], tfp1.df$end[i])
  tfp1.df$nmismatch[i] <- width(paste0(static,tfp1)) - sum(str_split(tfp1.df$match.seq[i], "")[[1]] == 
                                                             str_split(paste0(static,tfp1), "")[[1]])
  
  # y <- matchPattern(paste0(static,tfp1), DNAString(tfp1.df$seq[1]))
  # tfp1.df$match.seq[i] <- as.character(y[1])
  # tfp1.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
} }



#### Check for TFP2
x <- vmatchPattern(paste0(static,tfp2), reads, max.mismatch = 5)
tfp2.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp2.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp2.df$id])))
tfp2.df <- tfp2.df %>% left_join(y, by = "id")
tfp2.df$match.seq <- NA
tfp2.df$nmismatch <- NA
if(nrow(tfp2.df) > 0){
for(i in 1:nrow(tfp2.df)){ 
  tfp2.df$match.seq[i] <- str_sub(tfp2.df$seq[i], tfp2.df$start[i], tfp2.df$end[i])
  tfp2.df$nmismatch[i] <- width(paste0(static,tfp2)) - sum(str_split(tfp2.df$match.seq[i], "")[[1]] == 
                                str_split(paste0(static,tfp2), "")[[1]])

  # y <- matchPattern(paste0(static,tfp2), DNAString(tfp2.df$seq[1]))
  # tfp2.df$match.seq[i] <- as.character(y[1])
  # tfp2.df$nmismatch[i] <- nmismatch(paste0(static,tfp2), y)
  }}


 
#### Check for TFP3
x <- vmatchPattern(paste0(static,tfp3), reads, max.mismatch = 5)
tfp3.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp3.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp3.df$id])))
tfp3.df <- tfp3.df %>% left_join(y, by = "id")
tfp3.df$match.seq <- NA
tfp3.df$nmismatch <- NA
if(nrow(tfp3.df) > 0){
for(i in 1:nrow(tfp3.df)){ 
  tfp3.df$match.seq[i] <- str_sub(tfp3.df$seq[i], tfp3.df$start[i], tfp3.df$end[i])
  tfp3.df$nmismatch[i] <- width(paste0(static,tfp3)) - sum(str_split(tfp3.df$match.seq[i], "")[[1]] == 
                                                             str_split(paste0(static,tfp3), "")[[1]])
  
  # y <- matchPattern(paste0(static,tfp3), DNAString(tfp3.df$seq[1]))
  # tfp3.df$match.seq[i] <- as.character(y[1])
  # tfp3.df$nmismatch[i] <- nmismatch(paste0(static,tfp3), y)
}}




#### Check for TFP4
x <- vmatchPattern(paste0(static,tfp4), reads, max.mismatch = 5)
tfp4.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp4.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp4.df$id])))
tfp4.df <- tfp4.df %>% left_join(y, by = "id")
tfp4.df$match.seq <- NA
tfp4.df$nmismatch <- NA
if(nrow(tfp4.df) > 0){
for(i in 1:nrow(tfp4.df)){ 
  tfp4.df$match.seq[i] <- str_sub(tfp4.df$seq[i], tfp4.df$start[i], tfp4.df$end[i])
  tfp4.df$nmismatch[i] <- width(paste0(static,tfp4)) - sum(str_split(tfp4.df$match.seq[i], "")[[1]] == 
                                                             str_split(paste0(static,tfp4), "")[[1]])
  
  # y <- matchPattern(paste0(static,tfp4), DNAString(tfp4.df$seq[1]))
  # tfp4.df$match.seq[i] <- as.character(y[1])
  # tfp4.df$nmismatch[i] <- nmismatch(paste0(static,tfp4), y)
}}


# Clean up
rm(list = setdiff(ls(), c("tfp1.df", "tfp2.df", "tfp3.df", "tfp4.df")))


z1 <- min(width(tfp1.df$seq))
z2 <- min(width(tfp2.df$seq))
z3 <- min(width(tfp3.df$seq))
z4 <- min(width(tfp4.df$seq))


message(paste0("Pattern #1 min length match = ", z1))
message(paste0("Pattern #2 min length match = ", z2))
message(paste0("Pattern #3 min length match = ", z3))
message(paste0("Pattern #4 min length match = ", z4))


 #
 
########################################
# Graph results
########################################

rm(list = ls())
df <- read_tsv("junction_totals.v2.txt")
df$percent <- round(df$counts / df$total.reads.matched * 100,1)
df$category <- paste0(df$source, " | ", df$match.length*2, " bp query")
df <- df %>% group_by(category) %>% mutate(label_y = cumsum(percent) - 0.5 * percent) %>% ungroup()

ggplot(data = df, aes(x = category, y = percent, fill = junction)) +
  geom_bar(position = "stack", stat = "identity") + 
  #geom_text(aes(label = percent))
  geom_text(aes(label = round(counts,1)), position = position_stack(vjust = .5)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






#################################################################
# Now search the junctions that you think should be impossible
#################################################################


# Start clean
rm(list = ls())

# Import FASTQ from plasmidsaurus
#reads <- readDNAStringSet("../Sternberg_S2F_1/reads/raw_reads.fastq", format = "fastq")
#reads <- readDNAStringSet("SRR2127654.fastq", format = "fastq")
reads <- readDNAStringSet("SRR2135242.fastq", format = "fastq")


# Add the different fragments that you'll use to search
static <- "GACCGGGGCAACGAATGTCGCTGACGCTCGCACAAACCTC"
bad1 <- "GGTTTAGGAACATCAGCCATACTTAATGCGCGGTCCAACG"
bad2 <- "GGTTTGGTAGACAGCAATGGTTACGTGCCTGTGTCACTGG"
bad3 <- "GGTTTAGGAAGTAGCGCGACACGGGATGCTTACAGCTCGA"
bad4 <- "GACCTTTATTCACCCGCATCAGCAGTCATGGCCAGCTCTG"


# Try shortening the different fragments that you'll use to search
total.length <- 25 # length upstream and downstream of junction
static <- str_sub(static, -total.length, -1)
bad1 <- str_sub(bad1, 1, total.length)
bad2 <- str_sub(bad2, 1, total.length)
bad3 <- str_sub(bad3, 1, total.length)
bad4 <- str_sub(bad4, 1, total.length)

# Initiate function to take string and return rev. comp. as string
turn.me <- function(o){ return(as.character(reverseComplement(DNAString(o)))) }

#### Check for bad1
x <- vmatchPattern(paste0(turn.me(bad2),bad4), reads, max.mismatch = 5)
bad1.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% bad1.df$id],
            seq = unname(as.character(reads[names(reads) %in% bad1.df$id])))
bad1.df <- bad1.df %>% left_join(y, by = "id")
bad1.df$match.seq <- NA
bad1.df$nmismatch <- NA
if(nrow(bad1.df) > 0){
  for(i in 1:nrow(bad1.df)){ 
    bad1.df$match.seq[i] <- str_sub(bad1.df$seq[i], bad1.df$start[i], bad1.df$end[i])
    bad1.df$nmismatch[i] <- width(paste0(static,bad1)) - sum(str_split(bad1.df$match.seq[i], "")[[1]] == 
                                                               str_split(paste0(static,bad1), "")[[1]])
    
    # y <- matchPattern(paste0(static,bad1), DNAString(bad1.df$seq[1]))
    # bad1.df$match.seq[i] <- as.character(y[1])
    # bad1.df$nmismatch[i] <- nmismatch(paste0(static,bad1), y)
  } }



#### Check for bad2
x <- vmatchPattern(paste0(turn.me(bad3),bad4), reads, max.mismatch = 5)
bad2.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% bad2.df$id],
            seq = unname(as.character(reads[names(reads) %in% bad2.df$id])))
bad2.df <- bad2.df %>% left_join(y, by = "id")
bad2.df$match.seq <- NA
bad2.df$nmismatch <- NA
if(nrow(bad2.df) > 0){
  for(i in 1:nrow(bad2.df)){ 
    bad2.df$match.seq[i] <- str_sub(bad2.df$seq[i], bad2.df$start[i], bad2.df$end[i])
    bad2.df$nmismatch[i] <- width(paste0(static,bad2)) - sum(str_split(bad2.df$match.seq[i], "")[[1]] == 
                                                               str_split(paste0(static,bad2), "")[[1]])
    
    # y <- matchPattern(paste0(static,bad2), DNAString(bad2.df$seq[1]))
    # bad2.df$match.seq[i] <- as.character(y[1])
    # bad2.df$nmismatch[i] <- nmismatch(paste0(static,bad2), y)
  }}



#### Check for bad3
x <- vmatchPattern(paste0(turn.me(bad2),bad1), reads, max.mismatch = 5)
bad3.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% bad3.df$id],
            seq = unname(as.character(reads[names(reads) %in% bad3.df$id])))
bad3.df <- bad3.df %>% left_join(y, by = "id")
bad3.df$match.seq <- NA
bad3.df$nmismatch <- NA
if(nrow(bad3.df) > 0){
  for(i in 1:nrow(bad3.df)){ 
    bad3.df$match.seq[i] <- str_sub(bad3.df$seq[i], bad3.df$start[i], bad3.df$end[i])
    bad3.df$nmismatch[i] <- width(paste0(static,bad3)) - sum(str_split(bad3.df$match.seq[i], "")[[1]] == 
                                                               str_split(paste0(static,bad3), "")[[1]])
    
    # y <- matchPattern(paste0(static,bad3), DNAString(bad3.df$seq[1]))
    # bad3.df$match.seq[i] <- as.character(y[1])
    # bad3.df$nmismatch[i] <- nmismatch(paste0(static,bad3), y)
  }}




#### Check for bad4
x <- vmatchPattern(paste0(turn.me(bad3),bad1), reads, max.mismatch = 5)
bad4.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% bad4.df$id],
            seq = unname(as.character(reads[names(reads) %in% bad4.df$id])))
bad4.df <- bad4.df %>% left_join(y, by = "id")
bad4.df$match.seq <- NA
bad4.df$nmismatch <- NA
if(nrow(bad4.df) > 0){
  for(i in 1:nrow(bad4.df)){ 
    bad4.df$match.seq[i] <- str_sub(bad4.df$seq[i], bad4.df$start[i], bad4.df$end[i])
    bad4.df$nmismatch[i] <- width(paste0(static,bad4)) - sum(str_split(bad4.df$match.seq[i], "")[[1]] == 
                                                               str_split(paste0(static,bad4), "")[[1]])
    
    # y <- matchPattern(paste0(static,bad4), DNAString(bad4.df$seq[1]))
    # bad4.df$match.seq[i] <- as.character(y[1])
    # bad4.df$nmismatch[i] <- nmismatch(paste0(static,bad4), y)
  }}


# Clean up
rm(list = setdiff(ls(), c("bad1.df", "bad2.df", "bad3.df", "bad4.df")))


z1 <- min(width(bad1.df$seq))
z2 <- min(width(bad2.df$seq))
z3 <- min(width(bad3.df$seq))
z4 <- min(width(bad4.df$seq))


message(paste0("Pattern #1 min length match = ", z1))
message(paste0("Pattern #2 min length match = ", z2))
message(paste0("Pattern #3 min length match = ", z3))
message(paste0("Pattern #4 min length match = ", z4))


########################################
# Graph results
########################################

rm(list = ls())
df <- read_tsv("bad_junction_totals.txt")
df$percent <- round(df$counts / df$total.reads.matched * 100,1)
df$category <- paste0(df$source, " | ", df$match.length*2, " bp query")
df <- df %>% group_by(category) %>% mutate(label_y = cumsum(percent) - 0.5 * percent) %>% ungroup()

ggplot(data = df, aes(x = category, y = percent, fill = junction)) +
  geom_bar(position = "stack", stat = "identity") + 
  #geom_text(aes(label = percent))
  geom_text(aes(label = round(counts,1)), position = position_stack(vjust = .5)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))









#################################################################
# Okay, I missed 4 potentials at junction 2. Just redo this
#################################################################

# Start clean
rm(list = ls())

# Import FASTQ from plasmidsaurus
#reads <- readDNAStringSet("../Sternberg_S2F_1/reads/raw_reads.fastq", format = "fastq")
#reads <- readDNAStringSet("SRR2127654.fastq", format = "fastq")
reads <- readDNAStringSet("SRR2135242.fastq", format = "fastq")


# Add the different fragments that you'll use to search
tfp1 <- "GGTTTAGGAACATCAGCCATACTTAATGCGCGGTCCAACG"
tfp2 <- "GGTTTGGTAGACAGCAATGGTTACGTGCCTGTGTCACTGG"
tfp3 <- "GGTTTAGGAAGTAGCGCGACACGGGATGCTTACAGCTCGA"
tfp4 <- "GACCTTTATTCACCCGCATCAGCAGTCATGGCCAGCTCTG"

# Initiate function to take string and return rev. comp. as string
turn.me <- function(o){ return(as.character(reverseComplement(DNAString(o)))) }

# Store RCs
rc.tfp1 <- turn.me(tfp1)
rc.tfp2 <- turn.me(tfp2)
rc.tfp3 <- turn.me(tfp3)
rc.tfp4 <- turn.me(tfp4)


# Try shortening the different fragments that you'll use to search
total.length <- 40 # length upstream and downstream of junction
tfp1 <- str_sub(tfp1, 1, total.length)
tfp2 <- str_sub(tfp2, 1, total.length)
tfp3 <- str_sub(tfp3, 1, total.length)
tfp4 <- str_sub(tfp4, 1, total.length)
rc.tfp1 <- str_sub(rc.tfp1, 1, total.length)
rc.tfp2 <- str_sub(rc.tfp2, 1, total.length)
rc.tfp3 <- str_sub(rc.tfp3, 1, total.length)
rc.tfp4 <- str_sub(rc.tfp4, 1, total.length)


#### Check for Junction possibility 1
p <- paste0(rc.tfp4,tfp3)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp1.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp1.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp1.df$id])))
tfp1.df <- tfp1.df %>% left_join(y, by = "id")
tfp1.df$match.seq <- NA
tfp1.df$nmismatch <- NA
if(nrow(tfp1.df) > 0){
  for(i in 1:nrow(tfp1.df)){ 
    tfp1.df$match.seq[i] <- str_sub(tfp1.df$seq[i], tfp1.df$start[i], tfp1.df$end[i])
    tfp1.df$nmismatch[i] <- width(p) - sum(str_split(tfp1.df$match.seq[i], "")[[1]] == 
                                                               str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp1.df$seq[1]))
    # tfp1.df$match.seq[i] <- as.character(y[1])
    # tfp1.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }



#### Check for Junction possibility 2
p <- paste0(rc.tfp4,tfp2)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp2.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp2.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp2.df$id])))
tfp2.df <- tfp2.df %>% left_join(y, by = "id")
tfp2.df$match.seq <- NA
tfp2.df$nmismatch <- NA
if(nrow(tfp2.df) > 0){
  for(i in 1:nrow(tfp2.df)){ 
    tfp2.df$match.seq[i] <- str_sub(tfp2.df$seq[i], tfp2.df$start[i], tfp2.df$end[i])
    tfp2.df$nmismatch[i] <- width(p) - sum(str_split(tfp2.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp2.df$seq[1]))
    # tfp2.df$match.seq[i] <- as.character(y[1])
    # tfp2.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }



#### Check for Junction possibility 3
p <- paste0(rc.tfp2,tfp4)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp3.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp3.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp3.df$id])))
tfp3.df <- tfp3.df %>% left_join(y, by = "id")
tfp3.df$match.seq <- NA
tfp3.df$nmismatch <- NA
if(nrow(tfp3.df) > 0){
  for(i in 1:nrow(tfp3.df)){ 
    tfp3.df$match.seq[i] <- str_sub(tfp3.df$seq[i], tfp3.df$start[i], tfp3.df$end[i])
    tfp3.df$nmismatch[i] <- width(p) - sum(str_split(tfp3.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp3.df$seq[1]))
    # tfp3.df$match.seq[i] <- as.character(y[1])
    # tfp3.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }


#### Check for Junction possibility 4
p <- paste0(rc.tfp3,tfp4)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp4.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp4.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp4.df$id])))
tfp4.df <- tfp4.df %>% left_join(y, by = "id")
tfp4.df$match.seq <- NA
tfp4.df$nmismatch <- NA
if(nrow(tfp4.df) > 0){
  for(i in 1:nrow(tfp4.df)){ 
    tfp4.df$match.seq[i] <- str_sub(tfp4.df$seq[i], tfp4.df$start[i], tfp4.df$end[i])
    tfp4.df$nmismatch[i] <- width(p) - sum(str_split(tfp4.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp4.df$seq[1]))
    # tfp4.df$match.seq[i] <- as.character(y[1])
    # tfp4.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }


#### Check for Junction possibility 5
p <- paste0(rc.tfp1,tfp2)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp5.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp5.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp5.df$id])))
tfp5.df <- tfp5.df %>% left_join(y, by = "id")
tfp5.df$match.seq <- NA
tfp5.df$nmismatch <- NA
if(nrow(tfp5.df) > 0){
  for(i in 1:nrow(tfp5.df)){ 
    tfp5.df$match.seq[i] <- str_sub(tfp5.df$seq[i], tfp5.df$start[i], tfp5.df$end[i])
    tfp5.df$nmismatch[i] <- width(p) - sum(str_split(tfp5.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp5.df$seq[1]))
    # tfp5.df$match.seq[i] <- as.character(y[1])
    # tfp5.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }


#### Check for Junction possibility 6
p <- paste0(rc.tfp1,tfp3)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp6.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp6.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp6.df$id])))
tfp6.df <- tfp6.df %>% left_join(y, by = "id")
tfp6.df$match.seq <- NA
tfp6.df$nmismatch <- NA
if(nrow(tfp6.df) > 0){
  for(i in 1:nrow(tfp6.df)){ 
    tfp6.df$match.seq[i] <- str_sub(tfp6.df$seq[i], tfp6.df$start[i], tfp6.df$end[i])
    tfp6.df$nmismatch[i] <- width(p) - sum(str_split(tfp6.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp6.df$seq[1]))
    # tfp6.df$match.seq[i] <- as.character(y[1])
    # tfp6.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }

#### Check for Junction possibility 7
p <- paste0(rc.tfp2,tfp1)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp7.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp7.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp7.df$id])))
tfp7.df <- tfp7.df %>% left_join(y, by = "id")
tfp7.df$match.seq <- NA
tfp7.df$nmismatch <- NA
if(nrow(tfp7.df) > 0){
  for(i in 1:nrow(tfp7.df)){ 
    tfp7.df$match.seq[i] <- str_sub(tfp7.df$seq[i], tfp7.df$start[i], tfp7.df$end[i])
    tfp7.df$nmismatch[i] <- width(p) - sum(str_split(tfp7.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp7.df$seq[1]))
    # tfp7.df$match.seq[i] <- as.character(y[1])
    # tfp7.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }

#### Check for Junction possibility 8
p <- paste0(rc.tfp3,tfp1)
x <- vmatchPattern(p, reads, max.mismatch = 5)
tfp8.df <- tibble(id = names(x)[names(x) %in% row.names(data.frame(as.matrix(unlist(x))))],
                  start = unlist(startIndex(x)), end = unlist(endIndex(x)))
y <- tibble(id = names(reads)[names(reads) %in% tfp8.df$id],
            seq = unname(as.character(reads[names(reads) %in% tfp8.df$id])))
tfp8.df <- tfp8.df %>% left_join(y, by = "id")
tfp8.df$match.seq <- NA
tfp8.df$nmismatch <- NA
if(nrow(tfp8.df) > 0){
  for(i in 1:nrow(tfp8.df)){ 
    tfp8.df$match.seq[i] <- str_sub(tfp8.df$seq[i], tfp8.df$start[i], tfp8.df$end[i])
    tfp8.df$nmismatch[i] <- width(p) - sum(str_split(tfp8.df$match.seq[i], "")[[1]] == 
                                             str_split(p, "")[[1]])
    
    # y <- matchPattern(paste0(static,tfp1), DNAString(tfp8.df$seq[1]))
    # tfp8.df$match.seq[i] <- as.character(y[1])
    # tfp8.df$nmismatch[i] <- nmismatch(paste0(static,tfp1), y)
  } }


# Clean up
rm(list = setdiff(ls(), c("tfp1.df", "tfp2.df", "tfp3.df", "tfp4.df",
                          "tfp5.df", "tfp6.df", "tfp7.df", "tfp8.df")))


z1 <- min(width(tfp1.df$seq))
z2 <- min(width(tfp2.df$seq))
z3 <- min(width(tfp3.df$seq))
z4 <- min(width(tfp4.df$seq))
z5 <- min(width(tfp5.df$seq))
z6 <- min(width(tfp6.df$seq))
z7 <- min(width(tfp7.df$seq))
z8 <- min(width(tfp8.df$seq))



message(paste0("Pattern #1 min length match = ", z1))
message(paste0("Pattern #2 min length match = ", z2))
message(paste0("Pattern #3 min length match = ", z3))
message(paste0("Pattern #4 min length match = ", z4))
message(paste0("Pattern #5 min length match = ", z5))
message(paste0("Pattern #6 min length match = ", z6))
message(paste0("Pattern #7 min length match = ", z7))
message(paste0("Pattern #8 min length match = ", z8))


#

########################################
# Graph results
########################################

rm(list = ls())
df <- read_tsv("junction_totals.txt")
df$percent <- round(df$counts / df$total.reads.matched * 100,1)
df$category <- paste0(df$source, " | ", df$match.length*2, " bp query")
df <- df %>% group_by(category) %>% mutate(label_y = cumsum(percent) - 0.5 * percent) %>% ungroup()

ggplot(data = df, aes(x = category, y = percent, fill = junction)) +
  geom_bar(position = "stack", stat = "identity") + 
  #geom_text(aes(label = percent))
  geom_text(aes(label = round(counts,1)), position = position_stack(vjust = .5)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))






