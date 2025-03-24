# Volcano plots for 431
# 2/25/25

################################################
################ LOAD PACKAGES #################

library(ggplot2)
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
# library(knitr)
library(plotly)
# library(ggprism) # for add_pvalue()
# library(rstatix) # for adjust_pvalue
# library(ggpmisc) # https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
library(ggrepel)
# library(pheatmap)
# library(ggplotify) # To convert pheatmaps to ggplots
# library(corrplot)
# library(ggcorrplot)
# library(ggfortify) # To make pca plots with plotly



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
`W0.ComparedTo.Broth` <- read.delim("/Users/snork-maiden/Documents/Micro_grad_school/Sherman_Lab/R_projects/W0vsW2_GeneSetTesting/JOINED_BobAverages/MTb.MetaResults.W0_vs_Broth/W0.MTb.Meta.JOINED.txt")
`W2.ComparedTo.W0` <- read.delim("JOINED_BobAverages/MTb.MetaResults.W2_vs_W0/W2.MTb.Meta.JOINED.txt")
`W2.ComparedTo.Broth` <- read.delim("/Users/snork-maiden/Documents/Micro_grad_school/Sherman_Lab/R_projects/W0vsW2_GeneSetTesting/JOINED_BobAverages/MTb.MetaResults.W2_vs_Broth/W2.MTb.Meta.JOINED.txt")

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


