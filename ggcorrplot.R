# Correlation plot with ggcorrplot
# E. Lamont
# 4/3/25

source("Import_data.R") # to get my_tpm

# Log10 transform the data
my_tpm_Log10 <- my_tpm %>% 
  rename("W0_250754" = "S_250754",
         "W0_355466" = "S_355466",
         "W0_503557" = "S_503557",
         "W2_503937" = "S_503937",
         "W2_575533" = "S_575533_MtbrRNA",
         "W2_577208" = "S_577208") %>% 
  # mutate(Gene = rownames(my_tpm)) %>%
  mutate(across(where(is.numeric), ~ .x + 1)) %>% # Add 1 to all the values
  mutate(across(where(is.numeric), ~ log10(.x))) # Log transform the values


# Plot basics
my_plot_themes <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(10, 10, 10, 20),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_blank()
  )

# http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2

###########################################################
################# PEARSON LOG10 GGCORRPLOT ################

# Make the correlation
corr <- cor(my_tpm_Log10, method = "pearson")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_tpm_Log10)
# head(p.mat[, 1:4])

min(corr) # 0.2357629
my_min <- round(min(corr), 1)

# Plot pearson
ggcorrplot_PearsonLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             method = "square", 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  # scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Pearson Correlation Log10(TPM+1)", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_PearsonLog10

ggsave(ggcorrplot_PearsonLog10,
       file = "ggcorrplot_PearsonLog10_v3.pdf",
       path = "ggcorrplot_Figures",
       width = 7, height = 6, units = "in")


###########################################################
################ SPEARMAN LOG10 GGCORRPLOT ################
# Only going to do Pearson because the number of genes is so high CLT applies and parametric tests can be used

# Make the correlation
corr <- cor(my_tpm_Log10, method = "spearman")

# Compute a matrix of correlation p-values
p.mat <- cor_pmat(my_tpm_Log10)
# head(p.mat[, 1:4])

min(corr) # 0.2378325
my_min <- round(min(corr), 1)

# Plot
ggcorrplot_SpearmanLog10 <- corr %>% 
  ggcorrplot(hc.order = F, 
             lab = TRUE, lab_size = 4,
             type = c("lower"),
             outline.col = "white") + 
  scale_fill_gradient2(limit = c(my_min,1), low = "blue", high =  "red", mid = "white", midpoint = (((1-my_min)/2)+my_min)) + # Make sure to change based on the min!
  my_plot_themes + 
  scale_x_discrete(guide = guide_axis(angle = 45)) + 
  labs(title = "Spearman Correlation Log10(TPM+1)", 
       subtitle = NULL, 
       fill = "Correlation")
ggcorrplot_SpearmanLog10

# ggsave(ggcorrplot_SpearmanLog10,
#        file = "ggcorrplot_SpearmanLog10_v2.pdf",
#        path = "ggcorrplot_Figures",
#        width = 7, height = 6, units = "in")


###########################################################
################### TRYING OTHER THINGS ###################

# pairs(my_tpm_Log10)
# 
# # https://borisleroy.com/en/2013/06/09/correlation-plots-in-r/
# install.packages("Rarity")
# library(Rarity)

# Pearson
# pdf("ggcorrplot_Figures/rarity_PearsonLog10_v1.pdf", width = 10, height = 10)
# corPlot(my_tpm_Log10, method = "pearson", 
#         title = "Pearson Correlation Log10(TPM+1)") 
# dev.off()

# Spearman
# pdf("ggcorrplot_Figures/rarity_SpearmanLog10_v1.pdf", width = 10, height = 10)
# corPlot(my_tpm_Log10, method = "spearman", 
#         title = "Spearman Correlation Log10(TPM+1)") 
# dev.off()


