# Import Data with drug transcriptomes from prior work
# 6/5/25


###########################################################
#################### IMPORT RAW READS #####################

Sputum_w_Drug_rawReads <- read.csv("raw_data/INDIGO-transcriptomes_TestSputum_RawReads.csv")
Sputum_w_Drug_metadata <- read.csv("raw_data/INDIGO-transcriptomes_TestSputum_Metadata.csv")

# Remove rows that have NAs in them
Sputum_w_Drug_rawReads_2 <- Sputum_w_Drug_rawReads %>% 
  drop_na() %>% # Now 4027 rows
  column_to_rownames(var = "X") %>%
  as.matrix()


###########################################################
#################### BATCH CORRECTION #####################
# This basically all taken from Mark and I'm not sure what's happening
# https://academic.oup.com/nargab/article/2/3/lqaa078/5909519?login=true

batch <- Sputum_w_Drug_metadata$BatchNumber
Sputum_w_Drug_rawReads_BatchCorrected <- ComBat_seq(Sputum_w_Drug_rawReads_2, batch = batch)

# NOT WORKING!!
###########################################################
################# COVERT RAW READS TO TPM #################
## NOT WORKING@@ 


# Import Gene length info
load("MTb.MapSet.rda")
H37Rv_ExonMap <- mapSet$exonMap
H37Rv_GeneLengths <- H37Rv_ExonMap %>%
  mutate(GeneLength = END - POSITION) %>%
  select(GENE_ID, GeneLength)

rawReads_df <- test_df

# Make a rawCount -> TPM function
raw.to.tpm_func <- function(rawReads_df) {
  # TPM = (Read Count / Gene Length in kb) / sum(Read Count / Gene Length in kb for all genes) * 1e6

  # Ensure gene lengths are in the same order as raw_counts columns
  gene_lengths_ordered <- H37Rv_GeneLengths %>%
    filter(GENE_ID %in% rownames(rawReads_df)) %>%
    arrange(match(GENE_ID, rownames(rawReads_df)))  # Make sure order matches

  # Convert gene lengths to kilobases
  gene_lengths_kb <- (gene_lengths_ordered$GeneLength) / 1000

  # Calculate RPK (Reads Per Kilobase)
  rpk_df <- sweep(rawReads_df, 1, gene_lengths_kb, FUN = "/")

  # Calculate per-sample scaling factor (sum of RPKs/1e6). The per million scaling factor
  scaling_factors <- colSums(rpk_df)/1e6
  
  # Calculate TPM
  tpm_df <- sweep(rpk_df, 2, scaling_factors, FUN = "/")
}
test <- raw.to.tpm_func(Sputum_w_Drug_rawReads_2)
colSums(test) # check it's at 1 million
# These don't really match bob's...... Is it because so many of the genes got filtered out because they were not in the drug datasets?
test_df <- All_tpm %>% as.data.frame() %>% select("H37Ra_Broth_4", "S_355466")
test2 <- raw.to.tpm_func(test_df)
# No this still isn't right.....

###########################################################
#################### CONVERT TO CPM #######################

# Convert just the raw reads
Sputum_w_Drug_rawReads_cpm <- as.data.frame(cpm(Sputum_w_Drug_rawReads_2))

# Convert the batch corrected raw reads
Sputum_w_Drug_rawReads_BatchCorrected_cpm <- as.data.frame(cpm(Sputum_w_Drug_rawReads_BatchCorrected))















