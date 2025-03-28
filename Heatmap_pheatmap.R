# Try a heatmap with pheatmap
# E. Lamont 
# 3/25/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

# source("Import_data.R")
# I've just imported the DEG values now, not the tpm.... I guess just start with these, compare everything to broth
# No, this won't work because I need the replicates separated... need the tpm for all the samples...

# Start with my_tpm_W0vsBroth, also my_tpm which has all the unique sputum >1M and the broth

###########################################################
###################### PROCESS DATA #######################

# Add metadata for the GENEs!!
# Use Bob's MTb.TB.Phenotypes.AllGeneSets
### but what happens if one gene is in multiple sets?
my_tpm_W0vsBroth_2 <- my_tpm_W0vsBroth %>%
  rownames_to_column(var = "Gene") %>%
  mutate(GeneSet = sapply(Gene, function(g) {
    match_set <- names(MTb.TB.Phenotypes.AllGeneSets)[sapply(MTb.TB.Phenotypes.AllGeneSets, function(x) g %in% x)]
    if (length(match_set) > 0) match_set else "NA"
  }))


# Filter so there are fewer genes to deal with right now
my_tpm_W0vsBroth_filtered10 <- my_tpm_W0vsBroth %>%
  filter(if_all(where(is.numeric), ~ .x >= 10))  # Keep only rows where all columns are >=10
# Now there are only 4003 instead of 4499 rows.... doesn't help that much

# Lets filter for 100 just so it's easier for right now!
my_tpm_W0vsBroth_filtered100 <- my_tpm_W0vsBroth %>%
  filter(if_all(where(is.numeric), ~ .x >= 100))
# Now there are only 1341 rows! 

# Change to a matrix?
my_tpm_W0vsBroth_matrix <- as.matrix(my_tpm_W0vsBroth_filtered100)


###########################################################
######################## PHEATMAP #########################

# Not normalizing them because they are already in tpm, not sure if this is right....

pheatmap(my_tpm_W0vsBroth_matrix[1:10,], scale = "row")

# I need to scale by row to see anything...  not sure what is happening...

my_tpm_W0vsBroth_matrix[1:10,]


# Trying to subset based on gene set list

testing <- my_tpm_W0vsBroth %>% subset(rownames(my_tpm_W0vsBroth) %in% MTb.TB.Phenotypes.TopGeneSets$`human_sputum: top 50 genes`) # Guess this doesn't need to be a matrix

pheatmap(testing, scale = "row")

testing <- my_tpm %>% subset(rownames(my_tpm) %in% MTb.TB.Phenotypes.TopGeneSets$`human_sputum: top 50 genes`) # Guess this doesn't need to be a matrix

my_annotation_colors <- list(
  Week = c("0" = "#0072B2",  # Blue
           "2" = "#E66900",  # Orange
           "Broth" = "#999999")  # Grey
)

pheatmap(testing, 
         annotation_col = my_pipeSummary["Week"], 
         # annotation_row = gene_annot["Product"],
         annotation_colors = my_annotation_colors,
         scale = "row")


###########################################################
################### TESTING FOR SHINY #####################


allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["microaerophilic: top 25 genes"]]
my_data <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]])
p <- pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_colors = my_annotation_colors,
              scale = "row")
p
heatmap(as.matrix(my_data))
