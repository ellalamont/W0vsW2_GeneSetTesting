
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

# gene_annot <- read.delim("H37Rv.txt")
# row.names(gene_annot) <- gene_annot$Locus.Tag



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


###########################################################
############## IMPORT BOB's METAGENESETS DATA #############

# 6/18/25: Updated the txt files to include the analysis that was run after GSEA was fixed

`MetaGeneSets_W0vsBroth_UP` <- read.delim("MetaGeneSets_data/W0.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W0vsBroth_DOWN` <- read.delim("MetaGeneSets_data/W0.MTb.MetaGeneSets.DOWN.txt")
# Combine the UP and DOWN
# UP has a GSEA column, so removing that when binding
MetaGeneSets_W0vsBroth <- bind_rows(MetaGeneSets_W0vsBroth_UP, MetaGeneSets_W0vsBroth_DOWN) %>%
  select(intersect(names(MetaGeneSets_W0vsBroth_UP), names(MetaGeneSets_W0vsBroth_DOWN)))

# What about these other dataframes, are they more similar?
# `MetaGeneSets_BrothvsW0_UP` <- read.delim("MetaGeneSets_data/Broth.MTb.MetaGeneSets.UP.txt")
# `MetaGeneSets_BrothvsW0_DOWN` <- read.delim("MetaGeneSets_data/Broth.MTb.MetaGeneSets.DOWN.txt")

# MetaGeneSets_W0vsBroth_test <- bind_rows(MetaGeneSets_W0vsBroth_UP, MetaGeneSets_BrothvsW0_UP) %>%
#   select(intersect(names(MetaGeneSets_W0vsBroth_UP), names(MetaGeneSets_BrothvsW0_UP)))



# The W2 vs broth data
`MetaGeneSets_W2vsBroth_UP` <- read.delim("MetaGeneSets_data/W2.MTb.MetaGeneSets.UP.txt")
`MetaGeneSets_W2vsBroth_DOWN` <- read.delim("MetaGeneSets_data/W2.MTb.MetaGeneSets.DOWN.txt")
# Combine the UP and DOWN
# MetaGeneSets_W2vsBroth <- rbind(MetaGeneSets_W0vsBroth_UP, MetaGeneSets_W0vsBroth_DOWN)
# Error in rbind(deparse.level, ...) : 
#   numbers of columns of arguments do not match


