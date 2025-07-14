
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggtree)
library(DECIPHER)
library(rentrez)
library(rBLAST)

###################################################################################################################
# Use SRA genomic sequencing data to try to figure out junctions of Citrobacter prophage TFP genes

setwd("/Users/tannerwiegand/Documents/Projects/TnpB/dTnpB/Flagellin_manuscript/Figures/FigS4_assets/Citrobacter")
###################################################################################################################

rm(list = ls())

# Build BLAST database from Citrobacter reads
  # makeblastdb("SRR9220577.fasta")

# BLASTn search with suspected TFP-N region
  # bl <- blast(db = "SRR9220577.fasta")
  # 
  # query <- setNames(DNAStringSet(paste0("ATGGCACTCACACTATTGGCGGCCAACAACGCGCAAACGGTGCTGGCGGCAGGAATCAATGCAACTGCAACAA",
  #                                       "CTTTAACAGTTAACACAGGAACCGGTAATCTTTTCCCATCCCCAGTATCAGGGACCAGCTTTTTTAAACTGAC",
  #                                       "ACTTGTTGATTCTGCTACGGGCTCTCTCTCTGAAATTGTACATGTCACTTCCAGAACTGGTGACACGATGACC",
  #                                       "ATTGAGCGCGCCCAAGAGGGAACTACCGCGCGCATCTGGTCAGCAAATGACATTGCTGCAAACATGCTGACCG",
  #                                       "CAGGATCACTTCAACTCTACGCACAAAAAGATCAGTCTTTGTTGATTGCTAACAACCTTTCTGAGATAGCAAA",
  #                                       "TGCCGGACCGGATGCGGTTGCACAGACTCTCTCAAACC")),
  #                   "IR")
  # 
  # fmt <- paste("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
  #              "evalue", "length", "pident", "nident", "mismatch", "gaps", "qcov", "staxid")
  # 
  # cl <- predict(bl, query, custom_format = fmt)

# That didn't work

# Rename reads, since they have non-unique values and export
reads <- readDNAStringSet(c("SRR9220577.fasta",
                            "SRR9219931.fasta",
                            "SRR9220572.fasta"))
names(reads) <- paste0("read_", 1:length(reads))
# writeXStringSet(reads, "concat_reads.renamed.fna")
# makeblastdb("concat_reads.renamed.fna")

  # blastn -query TFP_N.fna -db concat_reads.renamed.fna -max_target_seqs 100000 -outfmt \
  # '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps' > \
  # TFP_N.concat_reads.blastout.E0.01.txt

colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
          "evalue", "length", "pident", "nident", "mismatch", "gaps", "qcov", "staxid")
blast <- read_tsv("TFP_N.concat_reads.blastout.E0.01.txt", col_names = colz)
  
# Add sequences
df <- tibble(id = names(reads), seq = unname(as.character(reads))) %>% distinct()
rm(reads)

blast <- blast %>% left_join(df %>% select(sseqid = id, seq))
rm(df)

### Extract the N-term region and any C-term regions
blast <- blast %>% mutate(upstream = NA, downstream = NA)

# Iterate through each hit to pull upstream and downstream sequences
for(i in 1:nrow(blast)){
  
  # If on positive strand
  if(blast$sstrand[i] == "plus"){
    
    # Extract upstream sequence, up to the start of the gene
    blast$upstream[i] <- str_sub(blast$seq[i], start = blast$sstart[i], end = blast$send[i])
    
    # Extract downstream sequence
    blast$downstream[i] <- str_sub(blast$seq[i], start = blast$send[i] + 1, end = -1)
    
  }
  
  # If on minus strand
  if(blast$sstrand[i] == "minus"){
    
    # Extract upstream sequence, up to the start of the gene
    blast$upstream[i] <- as.character(reverseComplement(DNAString(str_sub(blast$seq[i], start = blast$send[i], end = blast$sstart[i]))))
    
    # Extract downstream sequence
    blast$downstream[i] <- as.character(reverseComplement(DNAString(str_sub(blast$seq[i], start = 1, end = blast$send[i] - 1))))
    
  }
  
}

# Remove any hits that are shorter than 15 nt or where we didn't capture the downstream region
blast <- blast %>% filter(width(downstream) >= 15)

# Export an alignment of the downstream regions
blast <- blast %>% arrange(downstream, desc(width(downstream)))
writeXStringSet(setNames(DNAStringSet(blast$downstream), blast$sseqid), "Downstream_junction_regions.fna")

write_tsv(blast, "Final_BLASTout.tsv")




# Didn't end up using this strategy
# 
# # That's not working great. But I can manually find the IRs in a text editor
# 
# # Try string matching instead
# 
# # Search for hits
# hits.for <- unlist(vmatchPattern("CTCTCAAACC", reads, max.mismatch = 0)) # Pattern = predicted end of N-term
# 
# hits.rev <- unlist(vmatchPattern("CTCTCAAACC", reverseComplement(reads), max.mismatch = 0)) 
# 
# # Make a dataframe of hits
# hits.df <- rbind(df %>% filter(id %in% names(hits.for)) %>% mutate(strand = "+"),
#                  df.rev %>% filter(id %in% names(hits.rev)) %>% mutate(strand = "-"))
# 
# # Add start and end values
# x <- rbind(tibble(id = names(hits.for), start = start(hits.for), end = end(hits.for)),
#            tibble(id = names(hits.rev), start = start(hits.rev), end = end(hits.rev)))
# hits.df <- hits.df %>% left_join(x)
# 
# # Extract upstream and downstream sequences of the junction
# hits.df$upstream <- str_sub(hits.df$seq, start = hits.df$start, end = hits.df$end)
# hits.df$downstream <- str_sub(hits.df$seq, start = hits.df$end + 1, end = hits.df$end + 10)
# 
# downstream.count <- hits.df %>% filter(width(downstream) == 10) %>% count(downstream)
# upstream.count <- hits.df %>% filter(width(upstream) == 10) %>% count(upstream)
