# Volcano plots for 431
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
library(ggplotify) # To convert pheatmaps to ggplots
library(corrplot)
library(ggcorrplot)
library(ggfortify) # To make pca plots with plotly

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase")


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

# Altering these for the class! 

# load("MTb.KEGG.Pathways.rda")
# names(allGeneSets) <- gsub("<.*", "", names(allGeneSets))
# `Glycolysis / Gluconeogenesis ` <- allGeneSets$`Glycolysis / Gluconeogenesis ` # Pull this one out and will add it to the functional groups one
# 
# load("MTb.Tuberculist.FunctionalGroups.rda")
# # Add the glycolsis list
# allGeneSets <- append(allGeneSets, list(`Glycolysis / Gluconeogenesis ` = `Glycolysis / Gluconeogenesis `))

# Get list of all .rda files in the folder
rda_files <- list.files("GeneSet_Data", pattern = "\\.rda$", full.names = TRUE)

# Loop through each file and load it with a name based on the filename
for (file in rda_files) {
  file_name <- tools::file_path_sans_ext(basename(file))  # Extract filename without extension
  env <- new.env()  # Create a temporary environment to load the object
  load(file, envir = env)  # Load .rda file into this environment
  
  # Assign the loaded object to a new variable named after the file
  assign(file_name, env$allGeneSets)  
}

# Clean up
rm(env)  # Remove the temporary environment


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

my_tpm <- All_tpm %>% select(all_of(my_sample_names))
my_tpm_W0vsBroth <- All_tpm %>% select(all_of(W0vsBroth_sample_names))














