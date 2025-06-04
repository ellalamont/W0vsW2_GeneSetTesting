# Make a PCA plot with my samples of interest and Shuyi's drug transcriptomes
# 6/4/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R") # to get my_shuyi_tpm_filtered

# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=14),
        # legend.title = element_text(size = 14),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.text.x = element_text(angle = 0, size=14, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=14), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20)# ,
        # panel.background = element_rect(fill='transparent'),
        # plot.background = element_rect(fill='transparent', color=NA),
        # legend.background = element_rect(fill='transparent'),
        # legend.box.background = element_blank()
  )


###########################################################
############### MY_TPM WITH SHUYI DRUGS ALL ###############

# Transform the data
my_shuyi_tpm_filtered_t <- as.data.frame(t(my_shuyi_tpm_filtered))

# Remove columns that are all zero so the scale works for prcomp
my_shuyi_tpm_filtered_t2 <- my_shuyi_tpm_filtered_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_shuyi_tpm_filtered_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 29.5% of variance
summary_PCA[2,1] # PC2 explains 12.6% of variance
summary_PCA[3,1] # PC3 explains 10.0% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
# my_PCA_df <- merge(my_PCA_df, my_metadata, by = "SampleID", )

labels <- c(
  "S_250754" = "W0 sputum",
  "S_355466" = "W0 sputum",
  "S_503557" = "W0 sputum",
  "S_503937" = "W2 sputum",
  "S_575533_MtbrRNA" = "W2 sputum",
  "S_577208" = "W2 sputum",
  "H37Ra_Broth_4" = "Not captured broth",
  "H37Ra_Broth_5" = "Not captured broth",
  "H37Ra_Broth_6" = "Not captured broth",
  "5_EMBTRUE_D0_DrugRNAseq_U_R" = "EMB",
  "INH" = "Untreated",
  "PZA" = "Control",
  "RIF" = "Control"
)

# Add the labels I want by hand here, NEED TO DOUBLE CHECK THEY ARE CORRECT!
my_PCA_df <- my_PCA_df %>% 
  mutate(Labelling = c())


fig_PC1vsPC2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point() + 
  # geom_point(aes(fill = Labelling, shape = Labelling), size = 6, alpha = 0.8, stroke = 0.8) + 
  # scale_fill_manual(values=c(`W0 sputum` = "#0072B2", `W2 sputum` = "#E66900", `Not captured broth`= "#999999")) +  
  # scale_shape_manual(values=c(`W0 sputum` = 21, `W2 sputum` = 22, `Not captured broth`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_text(aes(label = SampleID)) + 
  labs(title = "PCA ",
       # subtitle = "All normal Depletion, no thresholds",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
# my_plot_themes_thumbnail
fig_PC1vsPC2
