
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(rBLAST)

###################################################################################################################
# Use SRA genomic sequencing data to try to figure out junctions of Citrobacter prophage TFP genes

setwd("path/to/Assess_Citrobacter_TFP")
###################################################################################################################

# Start by downloading the accompanying data from:
  # https://drive.google.com/drive/folders/1_UkuLmoi3Iri54BtBHwFQ73zhKYf8wC2?usp=drive_link

rm(list = ls())

# Manually downloaded reads from the NCBI SRA (as FASTA, no clipping/filtering)
  # SRR9220577.fasta
  # SRR9219931.fasta
  # SRR9220572.fasta

# Rename reads, since they have non-unique values and export
reads <- readDNAStringSet(c("SRR9220577.fasta",
                            "SRR9219931.fasta",
                            "SRR9220572.fasta"))
names(reads) <- paste0("read_", 1:length(reads))
writeXStringSet(reads, "concat_reads.renamed.fna")
makeblastdb("concat_reads.renamed.fna")

# BLASTn search of N-terminal fragment against reads

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

# Manually assessed and counted downstream junction regions

