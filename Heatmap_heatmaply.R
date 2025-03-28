# Try a heatmap with Heatmaply
# E. Lamont 
# 3/25/25

# https://cran.r-project.org/web/packages/heatmaply/vignettes/heatmaply.html

# source("Import_data.R")
# Start with my_tpm_W0vsBroth, also my_tpm which has all the unique sputum >1M and the broth


install.packages('heatmaply')
library("heatmaply")



###########################################################
######################## HEATMAPLY ########################

my_data <- my_tpm %>% subset(rownames(my_tpm) %in% allGeneSetList[["MTb.TB.Phenotypes.TopGeneSets"]][["human_sputum: top 25 genes"]])

heatmaply(my_data,
          scale = "row",
          col_side_colors = my_pipeSummary["Week"])
