# Try a heatmap with pheatmap
# E. Lamont 
# 3/25/25

# https://rpubs.com/tgjohnst/heatmaps_testing_1

# source("Import_data.R")
# I've just imported the DEG values now, not the tpm.... I guess just start with these, compare everything to broth
# No, this won't work because I need the replicates separated... need the tpm for all the samples...

# Start with my_tpm_W0vsBroth, also my_tpm which has all the unique sputum >1M and the broth

my_annotation_colors <- list(
  Week = c("0" = "#0072B2",  # Blue
           "2" = "#E66900",  # Orange
           "Broth" = "#999999")  # Grey
)

###########################################################
###################### PROCESS DATA #######################

# Add metadata for the GENEs!!
# Use Bob's MTb.TB.Phenotypes.AllGeneSets
### but what happens if one gene is in multiple sets?
# my_tpm_W0vsBroth_2 <- my_tpm_W0vsBroth %>%
#   rownames_to_column(var = "Gene") %>%
#   mutate(GeneSet = sapply(Gene, function(g) {
#     match_set <- names(MTb.TB.Phenotypes.AllGeneSets)[sapply(MTb.TB.Phenotypes.AllGeneSets, function(x) g %in% x)]
#     if (length(match_set) > 0) match_set else "NA"
#   }))


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


pheatmap(testing, 
         annotation_col = my_pipeSummary["Week"], 
         # annotation_row = gene_annot["Product"],
         annotation_colors = my_annotation_colors,
         scale = "row")


###########################################################
################# PHEATMAP FOR LANCE ######################
# 4/18/25: Lance want to check some of his genes of interest in my week 0 sputum samples

Lance_genes <- c("Rv3823c", "Rv1183", "Rv0206c", "Rv1886c", "Rv2942")
Lance_genes_2 <- c("Rv3820c", "Rv3821", "Rv3822", "Rv3823c", "Rv3824c", "Rv3825c", "Rv0295c", "Rv1182", "Rv3826")

Lance_subset <- my_tpm_W0vsBroth %>% subset(rownames(my_tpm_W0) %in% Lance_genes_2)

pheatmap(Lance_subset, 
         # annotation_row = Gene_Category, 
         annotation_col = my_pipeSummary["Week"],
         # annotation_colors = my_annotation_colors,
         scale = "row",
         display_numbers = T)

# write.csv(Lance_subset, "Lance_subset_2.csv")


###########################################################
######################## ALL DATA #########################

pheatmap(my_tpm , 
         # annotation_col = my_pipeSummary["Week"], 
         # annotation_row = gene_annot["Product"],
         annotation_colors = my_annotation_colors,
         scale = "row")

my_tpm_2 <- my_tpm[rowSums(my_tpm == 0) != ncol(my_tpm), ]

my_tpm_2_matrix <- my_tpm_2 %>% 
  # rename("W0_250754" = "S_250754",
  #        "W0_355466" = "S_355466",
  #        "W0_503557" = "S_503557",
  #        "W2_503937" = "S_503937",
  #        "W2_575533" = "S_575533_MtbrRNA",
  #        "W2_577208" = "S_577208") %>%
  as.matrix()

testing2 <- testing %>% 
  as.matrix()


pheatmap(testing2, 
         annotation_col = my_pipeSummary["Week"], 
         scale = "row",
         cutree_rows = 5,
         cutree_cols = 5)


###########################################################
############### ALL DATA WITH CLUSTERING ##################
# https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/

pheatmap(my_tpm_2_matrix, 
         # annotation_row = my_gene_col, 
         annotation_col = my_pipeSummary["Week"],
         scale = "row",
         cutree_rows = 5,
         cutree_cols = 5)

# Try to get good row annotations based on MTb functional group
# Start with MTb.TB.Phenotypes.AllGeneSets
# Convert to dataframe
Gene_Category <- do.call(rbind, lapply(names(Walter2015GeneSets), function(category) {
  data.frame(Gene = Walter2015GeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")


pheatmap(my_tpm_2_matrix, 
         # annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 8,
         cutree_cols = 5)





# Pull out what the clusters are: 
# use silent = TRUE to suppress the plot
my_heatmap <- pheatmap(my_tpm_2_matrix, 
                       annotation_row = Gene_Category, 
                       annotation_col = my_pipeSummary["Week"],
                       annotation_colors = my_annotation_colors,
                       scale = "row",
                       cutree_rows = 7,
                       cutree_cols = 5,
                       silent = TRUE)
# Extract the row clustering information
row_clusters <- cutree(my_heatmap$tree_row, k = 7)
# Convert to a data frame for easier handling
row_cluster_df <- data.frame(Gene = names(row_clusters), Cluster = row_clusters)
# View the first few rows
head(row_cluster_df)


###########################################################
############ W0 AND BROTH WITH CLUSTERING #################

my_tpm_3_matrix <- my_tpm_2 %>% select(-c("S_503937", "S_577208", "S_575533_MtbrRNA")) %>%
  as.matrix()

Gene_Category <- do.call(rbind, lapply(names(MTb.TB.Phenotypes.AllGeneSets), function(category) {
  data.frame(Gene = MTb.TB.Phenotypes.AllGeneSets[[category]], Category = category, stringsAsFactors = FALSE)
})) %>% 
  # filter(!Category %in% c("Cluster A", "Cluster B", "Cluster D", "Cluster E")) %>% 
  distinct(Gene, .keep_all = TRUE) %>% # Keep only the first gene occurance
  column_to_rownames(var = "Gene")

pheatmap(my_tpm_3_matrix, 
         annotation_row = Gene_Category, 
         fontsize_row = 1,
         annotation_col = my_pipeSummary["Week"],
         annotation_colors = my_annotation_colors,
         scale = "row",
         cutree_rows = 6,
         cutree_cols = 2)











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


selected_genes <- c("Rv0081", "Rv0494", "Rv2011c", "Rv1473A")
my_data <- my_tpm[rownames(my_tpm) %in% selected_genes, , drop = FALSE]
pheatmap(my_data, 
              annotation_col = my_pipeSummary["Week"], 
              annotation_row = gene_annot["Product"],  # Conditional annotation
              annotation_colors = my_annotation_colors,
              scale = "row", 
              fontsize = 18)



