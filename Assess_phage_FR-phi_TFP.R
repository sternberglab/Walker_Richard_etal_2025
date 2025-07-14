
library(tidyverse)
library(ggplot2)
library(Biostrings)

###################################################################################
# Check for reads that show evidence of tail fiber switching in BIDMC93 mutants (sSL3690)

setwd("path/to/Assess_phage_FR-ph_TFP")
###################################################################################

# Start by downloading the accompanying data from:
  # https://drive.google.com/drive/folders/1_UkuLmoi3Iri54BtBHwFQ73zhKYf8wC2?usp=drive_link

# Start clean
rm(list = ls())

# Make tibble with HTS file info
hts.data <- tibble(file = c("WTcmR-Prophage_HTS.fastq.gz",
                            "WT-Phage_HTS.fastq",
                            "delta_Gin_HTS.fastq.gz"), 
                   condition = c("WTcmR - Prophage (bacterial) sample",
                                 "WT - Phage sample",
                                 "∆gin"))

# Make tibble for results
df <- tibble()

# Iterate through each HTS dataset
for(g in 1:nrow(hts.data)){
  
  # Import FASTQ from plasmidsaurus
  reads <- readDNAStringSet(paste0(hts.data$file[g]), format = "fastq")
  
  # Add the different fragments that you'll use to search (These are hypothetical, based on manual genome viewing)
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
  
  
  
  
  #### Check for TFP4 [I expect this is impossible, since there's no inverted repeat to facilitate this configuration]
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
  rm(list = setdiff(ls(), c("tfp1.df", "tfp2.df", "tfp3.df", "tfp4.df", "hts.data", "df", 'g')))
  
  # Add to results tibble
  df <- rbind(df,
              tibble(species = c("N-C1", "N-C2", "N-C3", "N-C4"),
                     count = c(nrow(tfp1.df), nrow(tfp2.df), nrow(tfp3.df), nrow(tfp4.df)),
                     condition = hts.data$condition[g],
                     file = hts.data$file[g]
                     ))

}

# Export junction counts TSV
write_tsv(df, "Junction_totals.txt")

 
########################################
# Graph results
########################################

rm(list = ls())
df <- read_tsv("Junction_totals.txt")
df <- df %>% 
  group_by(condition) %>% 
  mutate(total.count = sum(count)) %>% 
  ungroup() %>% 
  mutate(percent = round(count / total.count * 100, 1)) %>% 
  group_by(condition) %>% 
  mutate(label_y = cumsum(percent) - 0.5 * percent) %>% 
  ungroup()

ggplot(data = df, aes(x = condition, y = percent, fill = species)) +
  geom_bar(position = "stack", stat = "identity") + 
  geom_text(aes(label = round(count,1)), position = position_stack(vjust = .5)) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


