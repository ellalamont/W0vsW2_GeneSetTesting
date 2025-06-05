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


