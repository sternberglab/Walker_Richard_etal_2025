
library(tidyverse)
library(ggplot2)
library(Biostrings)
library(ggtree)
library(DECIPHER)
library(rentrez)
library(rBLAST)

###################################################################################################################
# Annotate species on Matt's Gin/Hin tree

setwd("/Users/tannerwiegand/Documents/Projects/TnpB/dTnpB/Flagellin_manuscript/Figures/FigS4_assets")
###################################################################################################################

rm(list = ls())

# Import specices IDs (added by Matt)
meta <- read_tsv("Gin_tree_genus_ITOL.txt")

meta.count <- meta %>% count(genus)

# Make a list of only the genuses that have 10 or more representatives
these <- meta.count$genus[meta.count$n >= 10]

# Change other representatives to 0 count
meta$genus[!meta$genus %in% these] <- "Other"

# Also change the "candidatus" to other
meta$genus[meta$genus == "Candidatus"] <- "Other"

# Modify color of other category
meta$color[meta$genus == "Other"] <- "#A7A9AC"

# Fix IDs
meta$prot.acc <- gsub(" ", "_", meta$prot.acc)

# Export updated annotations
write_tsv(meta, "Gin_tree_genus_ITOL.v2.txt")

# Format genus counts for Illustrator graph
meta.count.bu <- meta.count
meta.count$genus[meta.count$genus == "Candidatus" | meta.count$n < 10] <- "Other"
meta.count <- meta.count %>% group_by(genus) %>% summarise(n = sum(n)) %>% arrange(desc(n))

cat(meta.count$genus, sep = "\t")
cat(meta.count$n, sep = "\t")

# Extract list of accessions and export
x <- gsub("\\'", "", meta$prot.acc)
write_tsv(tibble(id = x), "Gin_tree_accessions.txt", col_names = F)

# Pull IPG info and protein sequences
ipg <- tibble()
for(i in seq(1, length(x), by = 200)){
  
  # Post accession numbers
  these.acc <- x[i:(i+199)]
  search.token <- entrez_post(db = "protein", id = these.acc)
  Sys.sleep(1)
  
  # Fetch IPG info
  v <- entrez_fetch(db = "protein", rettype = "ipg", retmode = "text", web_history = search.token)
  Sys.sleep(1)
  
  # Convert to dataframe
  z <- read_tsv(v, show_col_types = FALSE)
  
  # Add to growing dataframe
  ipg <- rbind(ipg, z)
  
  # Fetch FASTA
  v <- entrez_fetch(db = "protein", rettype = "fasta", web_history = search.token)
  Sys.sleep(1)
  
  # Write to temporary file
  write(v, "temp.txt")
  
  # Import as xstringset
  z <- readAAStringSet("temp.txt")
  
  # Append to growing FASTA file and remove from memory/temp file
  writeXStringSet(z, "Gin_homologs_MW.faa", append = TRUE)
  rm(z)
  system("rm temp.txt")
  
  # Print status message
  message(paste0(round(((i+199)/length(x))*100, 2), "%..."))
}

# Run Hoodini
  # mamba activate flagcsnap
  # flagcsnap --input Gin_tree_accessions.txt --num-threads 6 --keep --api-key e6bf3c5ad8fe4dea2db927a8da490cdf3809

# Build BLAST database from neighborhoods
makeblastdb("neighborhoods.fasta")

  # tblastn -query Tn3_TPase.faa -db neighborhoods.fasta -evalue 0.01 -outfmt \
  # '6 qseqid sseqid qlen qstart qend sstart send sstrand evalue length pident nident mismatch gaps' > \
  # Tn3_TPase.Gin_neighborhoods.blastout.E0.01.txt

# Import BLAST results
colz <- c("qseqid", "sseqid", "qlen", "qstart", "qend", "sstart", "send", "sstrand",
          "evalue", "length", "pident", "nident", "mismatch", "gaps", "qcov", "staxid")
blast <- read_tsv("Tn3_TPase.Gin_neighborhoods.blastout.E0.01.txt", col_names = colz)

#######################################
# Combine Tn3 + TFP annotations
#######################################

rm(list = ls())

# Import ITOL associations (determined by MW and exported with the ITOL annotation spreadsheet)
colz <- c("id", "label", "color", "assoc")
tn3 <- read_tsv("Tn3_association.ITOL.MW.tsv", col_names = colz)
tfp <- read_tsv("TFP_association.ITOL.MW.tsv", col_names = colz)

# These should be listed in the same order, but double check this assumption
table(tn3$id == tfp$id) # 978 = TRUE. Looks good.

# Now see if we have any conflicting annotations, but doesn't look like it in ITOL
table(!is.na(tfp$assoc[tn3$assoc == "Tn3 associated"])) # 978 = FALSE. Looks good.

# Okay, now merge the annotations
itol <- tfp
itol$color[tn3$assoc == "Tn3 associated"] <- "#D05149"
itol$assoc[tn3$assoc == "Tn3 associated"] <- "Tn3 associated"

# Set white for the others
itol$color[is.na(itol$assoc)] <- "#FFFFFF"
itol$assoc[is.na(itol$assoc)] <- "None"

# Check annotations
table(itol$color)
table(itol$assoc)

# Looks good. Export
write_tsv(itol, "Tn3_TFP_assoc.concat.ITOL.MW.tsv")


#######################################
# Use Tn3 TPase HMM to annotate assoc
#######################################

rm(list = ls())

# Translate all ORFs in neighborhoods
  # getorf -sequence neighborhoods.fasta -outseq neighborhoods.ORFs.faa -table 1

# Used PHMMER on the Tn3 TPase sequence Matt sent me, then downloaded the HMM it hit to
  # Matt's TPase = Tn3_TPase_MW.faa
  # Resulting HMM = PF01526.hmm

# Use HMM to search ORFs
  # hmmsearch --tblout PF01526.neighborhoods_ORFs.hmmout.tbl --cut_ga PF01526.hmm neighborhoods.ORFs.faa

parse_hmmsearch_tbl <- function(file_path) {
  # Read all lines, skipping comments
  raw_lines <- readLines(file_path)
  data_lines <- raw_lines[!grepl("^#", raw_lines)]
  
  # Split lines by whitespace
  parsed_lines <- strsplit(trimws(data_lines), "\\s+")
  
  # Convert to a data frame
  tbl_df <- do.call(rbind, lapply(parsed_lines, function(x) {
    # Ensure at least 19 columns (some lines may lack the final 'FL=' field)
    length(x) <- max(19, length(x))
    x
  }))
  
  tbl_df <- as.data.frame(tbl_df, stringsAsFactors = FALSE)
  
  # Assign column names based on typical HMMER tbl output fields
  colnames(tbl_df) <- c(
    "target_name", "target_accession", 
    "query_name", "query_accession", 
    "full_E_value", "full_score", "full_bias",
    "best_1_domain_E_value", "best_1_domain_score", "best_1_domain_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc",
    "flags"
  )
  
  # Convert numeric columns
  numeric_cols <- c(
    "full_E_value", "full_score", "full_bias",
    "best_1_domain_E_value", "best_1_domain_score", "best_1_domain_bias",
    "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc"
  )
  
  tbl_df[numeric_cols] <- lapply(tbl_df[numeric_cols], function(x) suppressWarnings(as.numeric(x)))
  
  return(tbl_df)
}

# Import HMM results
hmm <- parse_hmmsearch_tbl("PF01526.neighborhoods_ORFs.hmmout.tbl")
hmm <- hmm[,1:19]

# Import ITOL annotations and the IPG info for these sequences
itol <- read_tsv("Tn3_TFP_assoc.concat.ITOL.MW.tsv")
ipg <- read_csv("records.csv")

# Add homolog info to HMMER output
hmm <- hmm %>% 
  mutate(nucleotide_id = word(target_name, start = 1, end = -4, sep = "_")) %>% 
  left_join(ipg %>% select(nucleotide_id, protein_id) %>% distinct())

# Check current annotations for sequences that had an HMM hit
table(itol$assoc[gsub("\\'", "", itol$id) %in% hmm$protein_id])
  # None        Tail fiber annotation        Tn3 associated 
  # 16                    12                    38 

# Okay, try changing just the ones that have no annotation currently
itol$assoc[gsub("\\'", "", itol$id) %in% hmm$protein_id & itol$assoc == "None"] <- "Tn3 associated"
itol$color[itol$assoc == "Tn3 associated"] <- "#D05149"

# Export as TSV
write_tsv(itol, "Tn3_TFP_assoc.concat.ITOL.MW.v2.tsv")
