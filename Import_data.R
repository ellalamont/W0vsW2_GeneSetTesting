
# 2/25/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(knitr)
library(plotly)
library(ggprism) # for add_pvalue()
library(rstatix) # for adjust_pvalue
library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
library(pheatmap)
# library(dendextend) # May need this for looking at pheatmap clustering
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly
library(edgeR) # for cpm
library(sva) # For ComBat_seq batch correction

# DuffyTools
library(devtools)
# install_github("robertdouglasmorrison/DuffyTools")
library(DuffyTools)
# install_github("robertdouglasmorrison/DuffyNGS")
# BiocManager::install("robertdouglasmorrison/DuffyTools")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biobase")


cbPalette_1 <- c("#999999", "#E69F00") # Gold and Grey
cbPalette_1.5 <- c("#E69F00", "#999999") # Gold and Grey
cbPalette_2 <- c( "#0072B2", "#999999") # Blue and Grey
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 <-  c("#bfbfbf", "#56B4E9")
cbPalette3 <-  c("#bfbfbf", "#E69F00")
cbPalette4 <- c("#56B4E9", "#009E73", "#F0E442")
c25 <- c(
  "dodgerblue2", "#E31A1C", "green4",
  "#6A3D9A","#FF7F00","black", "gold1",
  "skyblue2", "#FB9A99","palegreen2","#CAB2D6",
  "#FDBF6F","gray70", "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown"
)
c12 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black", "palegreen2", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4") 
c16 <- c("dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00", "black","gold1", "#FB9A99", "#CAB2D6", "palegreen2", "gray70", "maroon", "orchid1", "blue1", "darkturquoise", "darkorange4") 


# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default

###########################################################
################### IMPORT BOB's DE DATA ##################

# Don't know why these are only working with the full pathname....
`W0.ComparedTo.Broth` <- read.delim("JOINED_BobAverages/MTb.MetaResults.W0_vs_Broth/W0.MTb.Meta.JOINED.txt")
`W2.ComparedTo.W0` <- read.delim("JOINED_BobAverages/MTb.MetaResults.W2_vs_W0/W2.MTb.Meta.JOINED.txt")
`W2.ComparedTo.Broth` <- read.delim("JOINED_BobAverages/MTb.MetaResults.W2_vs_Broth/W2.MTb.Meta.JOINED.txt")

###########################################################
################ MAKE A LIST OF ALL DFs ###################
list_dfs <- list(`W0.ComparedTo.Broth`,
                 `W2.ComparedTo.W0`, 
                 `W2.ComparedTo.Broth`)

# Make a list of the names
df_names <- c("W0.ComparedTo.Broth",
              "W2.ComparedTo.W0", 
              "W2.ComparedTo.Broth")

# Give the df list the correct df names
names(list_dfs) <- df_names



###########################################################
############### ADD COLUMNS OF DE VALUES ##################

# Make a new list to hold dataframes with extra columns
list_dfs_2 <- list()

ordered_DE <- c("significant down", "not significant", "significant up")

# Add extra DE columns to each dataframe
for (i in 1:length(list_dfs)) {
  
  current_df <- list_dfs[[i]]
  current_df_name <- df_names[i]
  
  # Make the column pointing out which ones are differentially expressed
  current_df$DE <- ifelse(current_df$LOG2FOLD < -1 & current_df$AVG_PVALUE < 0.05, "significant down",
                          ifelse(current_df$LOG2FOLD > 1 & current_df$AVG_PVALUE < 0.05, "significant up", "not significant"))
  current_df$DE <- factor(current_df$DE, levels = ordered_DE)
  
  # Make the column with DE gene names for plotting on graph
  current_df$DE_labels <- ifelse(current_df$DE != "not significant", current_df$GENE_NAME, NA)
  
  list_dfs_2[[current_df_name]] <- current_df
  
}

###########################################################
###################### LOAD GENE SETS #####################

# To put them in a list of lists
rda_files <- list.files("GeneSet_Data", pattern = "\\.rda$", full.names = TRUE)
allGeneSetList <- list()

for(file in rda_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  env <- new.env()
  load(file, envir = env)  # loads allGeneSets into env
  allGeneSetList[[file_name]] <- env$allGeneSets  # store it in our list
  rm(env)
}

# Now update the names for each gene set in the list of lists:
allGeneSetList <- lapply(allGeneSetList, function(gset) {
  if (!is.null(gset)) {
    names(gset) <- gsub("<.*", "", names(gset))
  }
  return(gset)
})

###########################################################
############### IMPORT H37Rv GENE ANNOTATION ##############

gene_annot <- read.delim("H37Rv.txt")
row.names(gene_annot) <- gene_annot$Locus.Tag



###########################################################
############ IMPORT AND PROCESS ALL TPM VALUES ############
# NOT scaled

# ProbeTest5_tpm <- read.csv("ProbeTest5_Mtb.Expression.Gene.Data.SCALED.TPM.csv")
ProbeTest5_tpm <- read.csv("tpm_data/ProbeTest5_Mtb.Expression.Gene.Data.TPM_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_tpm <- read.csv("tpm_data/ProbeTest4_Mtb.Expression.Gene.Data.TPM.csv")
ProbeTest3_tpm <- read.csv("tpm_data/ProbeTest3_Mtb.Expression.Gene.Data.TPM.csv")

# Need to remove the undetermined which all share names
ProbeTest5_tpm$Undetermined_S0 <- NULL
ProbeTest4_tpm$Undetermined_S0 <- NULL
ProbeTest3_tpm$Undetermined_S0 <- NULL

# Merge the 3 documents
All_tpm <- merge(ProbeTest5_tpm, ProbeTest4_tpm, all = T)
All_tpm <- merge(All_tpm, ProbeTest3_tpm, all = T)

# Adjust the names so they are slightly shorter
names(All_tpm) <- gsub(x = names(All_tpm), pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

rownames(All_tpm) <- All_tpm[,1] # add the rownames


###########################################################
#################### SUBSET TPM VALUES ####################
# Just grab the sputum and broth samples that I want

# Unique Sputum: 
# W0 samples: "S_250754", "S_355466", "S_503557" 
# W2 samples: "S_349941_Probe_3D_25", "S_503937", "S_575533_MtbrRNA", "S_577208"
# W4 samples: "S_351946_Probe_4A_100", "S_687338_Probe_4A_100"

# Unique Sputum above 1M reads
# W0 samples: "S_250754", "S_355466", "S_503557" 
# W2 samples: "S_503937", "S_575533_MtbrRNA", "S_577208"

# Uncaptured Ra broth samples
# "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6"

my_sample_names <- c("S_250754", "S_355466", "S_503557", "S_503937", "S_575533_MtbrRNA", "S_577208", "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")
W0vsBroth_sample_names <- c("S_250754", "S_355466", "S_503557", "H37Ra_Broth_4", "H37Ra_Broth_5", "H37Ra_Broth_6")
W0_sample_names <- c("S_250754", "S_355466", "S_503557" )

my_tpm <- All_tpm %>% select(all_of(my_sample_names))
my_tpm_W0vsBroth <- All_tpm %>% select(all_of(W0vsBroth_sample_names))

my_tpm_W0 <- All_tpm %>% select(all_of(W0_sample_names))

###########################################################
############### IMPORT PIPELINE SUMMARY DATA ##############
# Importing the ProbeTests 3 and 4 and 5 to get all the sputum samples I have done

# ProbeTest5_pipeSummary <- read.csv("ProbeTest5_Pipeline.Summary.Details.csv")
# This has been edited to include more metadata!
ProbeTest5_pipeSummary <- read.csv("ProbeTest5_Pipeline.Summary.Details_moreTrim.csv") # This has the 3' end trimmed 40bp to increase the number of reads aligning
ProbeTest4_pipeSummary <- read.csv("ProbeTest4_Pipeline.Summary.Details.csv")
ProbeTest3_pipeSummary <- read.csv("ProbeTest3_Pipeline.Summary.Details.csv")

# Merge the 3 documents
All_pipeSummary <- merge(ProbeTest5_pipeSummary, ProbeTest4_pipeSummary, all = T)
All_pipeSummary <- merge(All_pipeSummary, ProbeTest3_pipeSummary, all = T)

All_pipeSummary$X <- NULL
All_pipeSummary$X.1 <- NULL

All_pipeSummary$Hyb_Time <- as.character(All_pipeSummary$Hyb_Time)
ordered_Hyb_Time <- c("4", "16")
All_pipeSummary$Hyb_Time <- factor(All_pipeSummary$Hyb_Time, levels = ordered_Hyb_Time)

All_pipeSummary$Week <- as.character(All_pipeSummary$Week)
ordered_Week <- c("0", "2", "4")
All_pipeSummary$Week <- factor(All_pipeSummary$Week, levels = ordered_Week)

All_pipeSummary$EukrRNADep <- as.character(All_pipeSummary$EukrRNADep)
ordered_EukrRNADep <- c("MtbrRNA", "DualrRNA")
All_pipeSummary$EukrRNADep <- factor(All_pipeSummary$EukrRNADep, levels = ordered_EukrRNADep)

# Remove the undetermined
All_pipeSummary <- All_pipeSummary %>% filter(SampleID != "Undetermined_S0")

# Remove the marmoset and the high low THP1 samples
All_pipeSummary <- All_pipeSummary %>% filter(!Sample_Type %in% c("Marmoset", "High_Low_THP1"))

All_pipeSummary$SampleID <- gsub(x = All_pipeSummary$SampleID, pattern = "_S.*", replacement = "") # This regular expression removes the _S and everything after it (I think...)

All_pipeSummary <- All_pipeSummary %>% mutate(Sputum_Number = str_extract(SampleID, "S_[0-9]+"))

# Just get the samples I am interested in
my_pipeSummary <- All_pipeSummary %>% filter(SampleID %in% (my_sample_names))
rownames(my_pipeSummary) <- my_pipeSummary$SampleID

# Change the NA to broth
my_pipeSummary$Week <- as.character(my_pipeSummary$Week)
my_pipeSummary$Week[is.na(my_pipeSummary$Week)] <- "Broth"
my_pipeSummary$Week <- as.factor(my_pipeSummary$Week)  # Convert back if needed

my_pipeSummary["Week"]


# Going in a different direction with importing in new sheet (6/5/25)
# ###########################################################
# ################# IMPORT SHUYI'S DRUG TXN #################
# 
# # This is all raw data, so will have to convert to TPM
# 
# ShuyiDrug_rawReads <- read.csv("raw_data/INDIGO-transcriptomes-all_v1.csv")
# ShuyiDrug_metadata <- read.csv("raw_data/INDIGO-metadata.csv")
# 
# ShuyiDrug_1 <- ShuyiDrug_rawReads %>% 
#   filter(Drug %in% c("EMB", "RIF", "PZA", "INH")) %>%
#   filter(Strain == "H37Rv")
# 
# 
# ###########################################################
# ################# COVERT RAW READS TO TPM #################
# 
# # Import Gene length info
# load("MTb.MapSet.rda")
# H37Rv_ExonMap <- mapSet$exonMap
# H37Rv_GeneLengths <- H37Rv_ExonMap %>%
#   mutate(GeneLength = END - POSITION) %>%
#   select(GENE_ID, GeneLength)
# 
# # Grab just the gene columns from the rawReads
# rawReadColumns_only <- ShuyiDrug_1 %>%
#   select(any_of(H37Rv_GeneLengths$GENE_ID))
# non_gene_metadata <- ShuyiDrug_1 %>%
#   select(-any_of(H37Rv_GeneLengths$GENE_ID))
# 
# # Make a rawCount -> TPM function
# raw.to.tpm_func <- function(rawReads_df) {
#   # TPM = (Read Count / Gene Length in kb) / sum(Read Count / Gene Length in kb for all genes) * 1e6
#   
#   # Ensure gene lengths are in the same order as raw_counts columns
#   gene_lengths_ordered <- H37Rv_GeneLengths %>%
#     filter(GENE_ID %in% colnames(rawReads_df)) %>%
#     arrange(match(GENE_ID, colnames(rawReads_df)))  # Make sure order matches
#   
#   # Convert gene lengths to kilobases
#   gene_lengths_kb <- (gene_lengths_ordered$GeneLength) / 1000
#   
#   # Calculate RPK (Reads Per Kilobase)
#   rpk_df <- sweep(rawReads_df, 2, gene_lengths_kb, FUN = "/")
#   
#   # Calculate per-sample scaling factor (sum of RPKs/1e6). The per million scaling factor
#   scaling_factors <- rowSums(rpk_df)/1e6
#   # Calculate TPM
#   tpm_df <- sweep(rpk_df, 1, scaling_factors, FUN = "/")
# }
# tpm_df <- raw.to.tpm_func(rawReadColumns_only)
# rowSums(tpm_df) # check it's at 1 million
# 
# # Add the metadata back to the tpm data
# ShuyiDrug_tpm <- cbind(non_gene_metadata, tpm_df)
# 
# # Average all the samples together
# ShuyiDrug_tpm_Average <- ShuyiDrug_tpm %>%
#   group_by(Drug) %>%
#   summarize(across(where(is.numeric), mean, na.rm = T)) %>%
#   select(-Batch) %>%
#   t() %>%
#   as.data.frame()
# colnames(ShuyiDrug_tpm_Average) <- ShuyiDrug_tpm_Average[1,] # Add column names
# ShuyiDrug_tpm_Average <- ShuyiDrug_tpm_Average[-1,] # Remove the old column name row
# ShuyiDrug_tpm_Average <- ShuyiDrug_tpm_Average %>%
#   mutate(across(everything(), as.numeric)) # Make everything numberic
# 
# # Do all the adjustments without the averages
# ShuyiDrug_tpm2 <- ShuyiDrug_tpm %>%
#   t() %>%
#   as.data.frame()
# colnames(ShuyiDrug_tpm2) <- ShuyiDrug_tpm2[1,] # Add column names
# ShuyiDrug_tpm2 <- ShuyiDrug_tpm2[-c(1:6),] # Remove the old column name row
# ShuyiDrug_tpm2 <- ShuyiDrug_tpm2 %>%
#   mutate(across(everything(), as.numeric)) # Make everything numberic
# 
# # ***********
# # Check the TPM is the same with my TPM samples
# ProbeTest5_rawReads <- read.csv("raw_data/ProbeTest5_Mtb.Expression.Gene.Data.ReadsM_moreTrim.csv")
# # transpose
# ProbeTest5_rawReads_t <- t(ProbeTest5_rawReads) %>% as.data.frame()
# colnames(ProbeTest5_rawReads_t) <- ProbeTest5_rawReads_t[1,]
# ProbeTest5_rawReads_t <- ProbeTest5_rawReads_t[-1,]
# ProbeTest5_rawReads_t <- ProbeTest5_rawReads_t %>%
#   mutate(across(everything(), as.numeric))
# ProbeTest5_rawReads_tpm_test <- raw.to.tpm_func(ProbeTest5_rawReads_t)
# ProbeTest5_rawReads_tpm_test_t <- t(ProbeTest5_rawReads_tpm_test)
# # Compare this to ProbeTest5_tpm... not identical, but issues may be due to rounding somewhere....
# # ***********
# 
# 
# ###########################################################
# ################ COMBINE MINE AND SHUYI TPM ###############
# # I'm using the TPM that Bob makes and the TPM I made from Shuyi's raw reads
# colSums(my_tpm) # Check this is the corrected TPM
# 
# my_shuyi_tpm <- merge(my_tpm, ShuyiDrug_tpm_Average, by = "row.names", all = TRUE)
# # Interesting that Shuyi's data is missing the MT and the rvn etc genes.....
# 
# # Keep only the genes that all samples have values for 
# my_shuyi_tpm_filtered <- merge(my_tpm, ShuyiDrug_tpm2, by = "row.names", all = F)
# rownames(my_shuyi_tpm_filtered) <- my_shuyi_tpm_filtered$Row.names
# my_shuyi_tpm_filtered <- my_shuyi_tpm_filtered[, -1]  # remove the now-redundant Row.names column
# 
# 
