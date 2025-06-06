# Make a PCA plot with my samples of interest and Shuyi's drug transcriptomes
# 6/4/25

# Look into ggbiplot for more PCA stuff??
# https://cran.r-project.org/web/packages/ggbiplot/readme/README.html

# http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

# Two options in base R, prcomp() and princomp()
# prcomp() is preferred according to the website above

source("Import_data.R") 
source("Import_data_w_Drugs.R") # To get Sputum_w_Drug_rawReads_BatchCorrected_cpm and Sputum_w_Drug_metadata and Sputum_w_Drug_rawReads_cpm

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
################### BATCH CORRECTED CPM ###################

# Transform the data
my_cpm_t <- as.data.frame(t(Sputum_w_Drug_rawReads_BatchCorrected_cpm))

# Remove columns that are all zero so the scale works for prcomp
my_cpm_t2 <- my_cpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_cpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 26.2% of variance
summary_PCA[2,1] # PC2 explains 13.0% of variance
summary_PCA[3,1] # PC3 explains 10.8% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, Sputum_w_Drug_metadata, by = "SampleID")


fig_PC1vsPC2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Label, shape = Label2), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Week 2 sputum` = "#E66900", `Log broth`= "#999999", `RIF` = "#D32F2F", `PZA` = "#008080", `INH` = "#808000", `EMB` = "#B39DDB")) +  
  scale_shape_manual(values=c(`Sputum` = 21, `Log broth` = 22, `Drug`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  geom_text(aes(label = SampleID)) + 
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
    shape = guide_legend(order = 1)) +
  labs(fill = "Sample Type", 
    shape = "Sample Type",
    title = "PCA: Batch Correct CPM, Sputum vs Drug transcriptomes",
    x = paste0("PC1: ", summary_PCA[1,1], "%"),
    y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2

ggsave(fig_PC1vsPC2,
       file = "PCA_SputumVsDrug_BatchCorrect_cpm.pdf",
       path = "PCA_Figures",
       width = 9, height = 6, units = "in")



###########################################################
######################## CPM ONLY #########################

# Transform the data
my_cpm_t <- as.data.frame(t(Sputum_w_Drug_rawReads_cpm))

# Remove columns that are all zero so the scale works for prcomp
my_cpm_t2 <- my_cpm_t %>% select_if(colSums(.) != 0)

# Make the actual PCA
my_PCA <- prcomp(my_cpm_t2, scale = TRUE)

# See the % Variance explained
summary(my_PCA)
summary_PCA <- format(round(as.data.frame(summary(my_PCA)[["importance"]]['Proportion of Variance',]) * 100, digits = 1), nsmall = 1) # format and round used to control the digits after the decimal place
summary_PCA[1,1] # PC1 explains 28.9% of variance
summary_PCA[2,1] # PC2 explains 11.4% of variance
summary_PCA[3,1] # PC3 explains 10.1% of variance

# MAKE PCA PLOT with GGPLOT 
my_PCA_df <- as.data.frame(my_PCA$x[, 1:3]) # Extract the first 3 PCs
my_PCA_df <- data.frame(SampleID = row.names(my_PCA_df), my_PCA_df)
my_PCA_df <- merge(my_PCA_df, Sputum_w_Drug_metadata, by = "SampleID")


fig_PC1vsPC2 <- my_PCA_df %>% 
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = Label, shape = Label2), size = 6, alpha = 0.8, stroke = 0.8) + 
  scale_fill_manual(values=c(`Week 0 sputum` = "#0072B2", `Week 2 sputum` = "#E66900", `Log broth`= "#999999", `RIF` = "#D32F2F", `PZA` = "#008080", `INH` = "#808000", `EMB` = "#B39DDB")) +  
  scale_shape_manual(values=c(`Sputum` = 21, `Log broth` = 22, `Drug`= 23)) + 
  # geom_text_repel(aes(label = Week), size= 2.5, box.padding = 0.4, segment.color = NA, max.overlaps = Inf) + 
  # geom_text(aes(label = SampleID)) + 
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21)),
         shape = guide_legend(order = 1)) +
  labs(fill = "Sample Type", 
       shape = "Sample Type",
       title = "PCA: CPM, Sputum vs Drug transcriptomes (No batch correction)",
       x = paste0("PC1: ", summary_PCA[1,1], "%"),
       y = paste0("PC2: ", summary_PCA[2,1], "%")) +
  my_plot_themes
fig_PC1vsPC2

ggsave(fig_PC1vsPC2,
       file = "PCA_SputumVsDrug_cpm.pdf",
       path = "PCA_Figures",
       width = 9, height = 6, units = "in")






