# Looking at the MetaGeneSets from Bob's pipeline and subsetting after everything has been run
# E. Lamont
# 6/13/25

# Here the gene set enrichment analysis has already been done in Bob's meta way and I am just visualizing the result


# I think there are more significant things when everything is run together but I am not sure why


source("Import_data.R") # To get MetaGeneSets_W0vsBroth

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "right",legend.text=element_text(size=9),
        legend.title = element_text(size = 10),
        # legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_text(size=10), 
        axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=10), 
        plot.margin = margin(10, 10, 10, 20)# ,
  )



###########################################################
################ iMODULONS: SUBSET THE DATA ###############

# Why are the UP and DOWN just a little bit different?
# I don't really understand this... I'll edit so the up ones are only the positive LOG2FOLD and the DOWN are only the negative LOG2FOLD. Is this right?
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_iModulons <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  

MetaGeneSets_W0vsBroth_DOWN_iModulons <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD <= 0) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  

# There are 72 iModulons accounted for now, so lost some somewhere.....

# Combine the two 
MetaGeneSets_W0vsBroth_ALL_iModulons <- rbind(MetaGeneSets_W0vsBroth_UP_iModulons %>% filter(LOG2FOLD >= 0), MetaGeneSets_W0vsBroth_DOWN_iModulons %>% filter(LOG2FOLD <= 0))

###########################################################
############## iMODULONS: MAKE A COLUMN PLOT ##############

# The UP file
ColumnPlot <- MetaGeneSets_W0vsBroth_UP_iModulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="LOG2FOLD",
       title="iModulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.UP.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.UP.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 18, units = "in")

# The DOWN file
ColumnPlot <- MetaGeneSets_W0vsBroth_DOWN_iModulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="LOG2FOLD",
       title="iModulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.DOWN.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.DOWN.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 18, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-3.5, 3.5), expand = c(0, 0)) +
  labs(x="iModulon", y="LOG2FOLD",
       title="iModulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 18, units = "in")

# ALL just significant
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% 
  filter(AVG_PVALUE < 0.05) %>%
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-3.5, 3.5), expand = c(0, 0)) +
  labs(x="iModulon", y="LOG2FOLD",
       title="iModulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=14))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_SignificantOnly.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 9, units = "in")


###########################################################
################ REGULONS: SUBSET THE DATA ################

# Why are the UP and DOWN just a little bit different?
# I don't really understand this... I'll edit so the up ones are only the positive LOG2FOLD and the DOWN are only the negative LOG2FOLD. Is this right?
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_Regulons <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  

MetaGeneSets_W0vsBroth_DOWN_Regulons <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim())  

# Combine the two 
MetaGeneSets_W0vsBroth_ALL_Regulons <- rbind(MetaGeneSets_W0vsBroth_UP_Regulons %>% filter(LOG2FOLD >= 0), 
                                              MetaGeneSets_W0vsBroth_DOWN_Regulons %>% filter(LOG2FOLD <= 0))


###########################################################
############## REGULONS: MAKE A COLUMN PLOT ###############

# The UP file
ColumnPlot <- MetaGeneSets_W0vsBroth_UP_Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Regulons", y="LOG2FOLD",
       title="Regulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.UP.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Regulons_W0.MTb.MetaGeneSets.UP.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")

# The DOWN file
ColumnPlot <- MetaGeneSets_W0vsBroth_DOWN_Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Regulons", y="LOG2FOLD",
       title="Regulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.DOWN.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Regulons_W0.MTb.MetaGeneSets.DOWN.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
  labs(x="Regulons", y="LOG2FOLD",
       title="Regulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Regulons_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")


# ALL just significant
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_Regulons %>% 
  filter(AVG_PVALUE < 0.05) %>%
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
  labs(x="Regulons", y="LOG2FOLD",
       title="Regulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=14))
ColumnPlot
ggsave(ColumnPlot,
       file = "Regulons_W0.MTb.MetaGeneSets.ALL_SignificantOnly.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 10, height = 10, units = "in")



###########################################################
############# TFOE.REGULONS: SUBSET THE DATA ##############

# Why are the UP and DOWN just a little bit different?
# I don't really understand this...
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_TFOE.Regulons <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "TFOE.Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim()) %>% 
  mutate(PathName = str_remove(PathName, "^TFOE\\.Regulons: ") %>% str_trim()) # Removes the beginning with TFOE.Regulons

MetaGeneSets_W0vsBroth_DOWN_TFOE.Regulons <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "TFOE.Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim()) %>% 
  mutate(PathName = str_remove(PathName, "^TFOE\\.Regulons: ") %>% str_trim()) # Removes the beginning with TFOE.Regulons

# Combine the two 
MetaGeneSets_W0vsBroth_ALL_TFOE.Regulons <- rbind(MetaGeneSets_W0vsBroth_UP_TFOE.Regulons %>% filter(LOG2FOLD >= 0), 
                                                  MetaGeneSets_W0vsBroth_DOWN_TFOE.Regulons %>% filter(LOG2FOLD <= 0))


###########################################################
########### TFOE.REGULONS: MAKE A COLUMN PLOT #############

# The UP file
ColumnPlot <- MetaGeneSets_W0vsBroth_UP_TFOE.Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -4.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
  labs(x="TFOE.Regulons", y="LOG2FOLD",
       title="TFOE.Regulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.UP.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "TFOE.Regulons_W0.MTb.MetaGeneSets.UP.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")

# The DOWN file
ColumnPlot <- MetaGeneSets_W0vsBroth_DOWN_TFOE.Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -4.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
  labs(x="TFOE.Regulons", y="LOG2FOLD",
       title="TFOE.Regulons in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.DOWN.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "TFOE.Regulons_W0.MTb.MetaGeneSets.DOWN.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_TFOE.Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -4.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
  labs(x="TFOE.Regulons", y="LOG2FOLD",
       title="TFOE.Regulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "TFOE.Regulons_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 30, units = "in")


# ALL just significant
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_TFOE.Regulons %>% 
  filter(AVG_PVALUE < 0.05) %>%
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -4.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4.5, 4.5), expand = c(0, 0)) +
  labs(x="TFOE.Regulons", y="LOG2FOLD",
       title="TFOE.Regulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=14))
ColumnPlot
ggsave(ColumnPlot,
       file = "TFOE.Regulons_W0.MTb.MetaGeneSets.ALL_SignificantOnly.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 10, height = 10, units = "in")


###########################################################
################# Tuberculist.GO.Ontology #################

# Why are the UP and DOWN just a little bit different?
# I don't really understand this...
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_Tuberculist.GO.ONTOLOGY <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "Tuberculist.GO.Ontology")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim()) 

MetaGeneSets_W0vsBroth_DOWN_Tuberculist.GO.ONTOLOGY <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "Tuberculist.GO.Ontology")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = PathName %>%
           str_replace("<.*", "") %>%        # remove anything after <
           str_remove_all("&nbsp;") %>%      # remove all &nbsp;
           str_trim()) 

# Combine the two 
MetaGeneSets_W0vsBroth_ALL_Tuberculist.GO.Ontology <- rbind(MetaGeneSets_W0vsBroth_UP_Tuberculist.GO.ONTOLOGY %>% filter(LOG2FOLD >= 0), 
                                                            MetaGeneSets_W0vsBroth_DOWN_Tuberculist.GO.ONTOLOGY %>% filter(LOG2FOLD <= 0))


# Column Plot The UP file
ColumnPlot <- MetaGeneSets_W0vsBroth_UP_Tuberculist.GO.ONTOLOGY %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.ONTOLOGY", y="LOG2FOLD",
       title="Tuberculist.GO.ONTOLOGY in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.UP.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Tuberculist.GO.ONTOLOGY_W0.MTb.MetaGeneSets.UP.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 49.9, units = "in")

# The DOWN file
ColumnPlot <- MetaGeneSets_W0vsBroth_DOWN_Tuberculist.GO.ONTOLOGY %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.ONTOLOGY", y="LOG2FOLD",
       title="Tuberculist.GO.ONTOLOGY in W0 Sputum vs Log broth using W0.MTb.MetaGeneSets.DOWN.txt",
       subtitle = NULL) + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Tuberculist.GO.ONTOLOGY_W0.MTb.MetaGeneSets.DOWN.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 49.9, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_Tuberculist.GO.Ontology %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.ONTOLOGY", y="LOG2FOLD",
       title="Tuberculist.GO.ONTOLOGY in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Tuberculist.GO.ONTOLOGY_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 49.9, units = "in")


# ALL just significant
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_Tuberculist.GO.Ontology %>% 
  filter(AVG_PVALUE < 0.05) %>%
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  geom_text(aes(y = -3.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-4, 4), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.ONTOLOGY", y="LOG2FOLD",
       title="Tuberculist.GO.ONTOLOGY in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=14))
ColumnPlot
ggsave(ColumnPlot,
       file = "Tuberculist.GO.ONTOLOGY_W0.MTb.MetaGeneSets.ALL_SignificantOnly.pdf",
       path = "GSEA_Figures/MetaGeneSets",
       width = 18, height = 7, units = "in")

###########################################################
############## EXTRACT ALL SIGNIFICANT VALUES #############

merged <- rbind(MetaGeneSets_W0vsBroth_ALL_iModulons, MetaGeneSets_W0vsBroth_ALL_Regulons, MetaGeneSets_W0vsBroth_ALL_TFOE.Regulons, MetaGeneSets_W0vsBroth_ALL_Tuberculist.GO.Ontology) %>% 
  filter(AVG_PVALUE < 0.05)

write.csv(merged, file = "GSEA_Figures/MetaGeneSets/MetaGeneSets_merged_SignficantlyOnly.csv")




