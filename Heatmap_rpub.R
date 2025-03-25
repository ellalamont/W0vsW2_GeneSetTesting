# Try a heatmap with gplot rpub
# E. Lamont 
# 3/25/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

# source("Import_data.R")
# I've just imported the DEG values now, not the tpm.... I guess just start with these, compare everything to broth
# No, this won't work because I need the replicates separated... need the tpm for all the samples...

# Start with my_tpm_W0vsBroth

###########################################################
###################### PROCESS DATA #######################

# Add metadata for the GENEs!!
# Use Bob's MTb.TB.Phenotypes.AllGeneSets
df <- my_tpm_W0vsBroth %>%
  rownames_to_column(var = "Gene") %>%
  mutate(GeneSet = sapply(Gene, function(g) {
    match_set <- names(MTb.TB.Phenotypes.AllGeneSets)[sapply(MTb.TB.Phenotypes.AllGeneSets, function(x) g %in% x)]
    if (length(match_set) > 0) match_set else "NA"
  }))


# Filter so there are fewer genes to deal with right now
my_tpm_W0vsBroth_filtered10 <- my_tpm_W0vsBroth %>%
  filter(if_all(everything(), ~ .x >= 10))  # Keep only rows where all columns are >=10
# Now there are only 4003 instead of 4499 rows.... doesn't help that much

# Lets filter for 100 just so it's easier for right now!
my_tpm_W0vsBroth_filtered100 <- my_tpm_W0vsBroth %>%
  filter(if_all(everything(), ~ .x >= 100))
# Now there are only 1341 rows! 

# Change to a matrix?
my_tpm_W0vsBroth_matrix <- as.matrix(my_tpm_W0vsBroth_filtered10)


###########################################################
####################### BASE HEATMAP ######################

heatmap(my_tpm_W0vsBroth_matrix)

# Don't cluster
heatmap(my_tpm_W0vsBroth_matrix,
        Rowv=NA, Colv=NA)

# Change the color
heatmap(my_tpm_W0vsBroth_matrix,
        col=rev(brewer.pal(9,"RdBu")))



###########################################################
###################### Gplot HEATMAP ######################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

# Default unscaled heatmap
heatmap.2(my_tpm_W0vsBroth_matrix, col=rev(brewer.pal(9,"RdBu")))

# Scaled within rows, no column clustering
heatmap.2(my_tpm_W0vsBroth_matrix, col=rev(brewer.pal(9,"RdBu")), scale="row", Colv = "NA")


if (!require("NMF")) {
  install.packages("NMF", dependencies = TRUE)
  library(NMF)
}
# Non-scaled
aheatmap(my_tpm_W0vsBroth_matrix, color = "-RdBu", annColors = "Set2", annRow=geneExp$type)







