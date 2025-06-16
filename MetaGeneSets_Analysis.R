# Looking at the MetaGeneSets from Bob's pipeline and subsetting after everything has been run
# E. Lamont
# 6/13/25


# I think there are more significant things when everything is run together but I am not sure why


source("Import_data.R") # To get MetaGeneSets_W0vsBroth


###########################################################
################ iMODULONS: SUBSET THE DATA ###############

# Why are the UP and DOWN just a little bit different?
# I don't really understand this... I'll edit so the up ones are only the positive LOG2FOLD and the DOWN are only the negative LOG2FOLD. Is this right?
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_iModulons <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD >= 0) %>% 
  mutate(PathName = str_replace(PathName, "<.*", "") %>% str_trim())

MetaGeneSets_W0vsBroth_DOWN_iModulons <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "iModulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  # filter(LOG2FOLD <= 0) %>%
  mutate(PathName = str_replace(PathName, "<.*", "") %>% str_trim())

# There are 72 iModulons accounted for now, so lost some somewhere.....

# Combine the two 
MetaGeneSets_W0vsBroth_ALL_iModulons <- rbind(MetaGeneSets_W0vsBroth_UP_iModulons %>% filter(LOG2FOLD >= 0), 
                                              MetaGeneSets_W0vsBroth_DOWN_iModulons %>% filter(LOG2FOLD <= 0))

# Shorten the PathName
# MetaGeneSets_W0vsBroth_ALL_iModulons$PathName <- str_replace(MetaGeneSets_W0vsBroth_ALL_iModulons$PathName, "^(([^:]*):([^:]*)):.*", "\\1")
# MetaGeneSets_W0vsBroth_ALL_iModulons$PathName <- sub("<.*", "", MetaGeneSets_W0vsBroth_ALL_iModulons$PathName)

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
       path = "GSEA_Figures",
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
       path = "GSEA_Figures",
       width = 18, height = 18, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="LOG2FOLD",
       title="iModulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures",
       width = 18, height = 18, units = "in")

###########################################################
################ REGULONS: SUBSET THE DATA ################

# Why are the UP and DOWN just a little bit different?
# I don't really understand this... I'll edit so the up ones are only the positive LOG2FOLD and the DOWN are only the negative LOG2FOLD. Is this right?
# Or lets just plot them separately, not LOG2FOLD filtering

MetaGeneSets_W0vsBroth_UP_Regulons <- MetaGeneSets_W0vsBroth_UP %>% 
  filter(str_detect(PathName, "Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = str_replace(PathName, "<.*", "") %>% str_trim())

MetaGeneSets_W0vsBroth_DOWN_Regulons <- MetaGeneSets_W0vsBroth_DOWN %>% 
  filter(str_detect(PathName, "Regulons")) %>% 
  filter(!str_detect(PathName, "ISB.Corems")) %>%
  filter(!str_detect(PathName, "TFOE.Regulons")) %>%
  mutate(PathName = str_replace(PathName, "<.*", "") %>% str_trim())

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
       path = "GSEA_Figures",
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
       path = "GSEA_Figures",
       width = 18, height = 30, units = "in")

# The two files combined 
ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_Regulons %>% 
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  # geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  # scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Regulons", y="LOG2FOLD",
       title="Regulons in W0 Sputum vs Log broth using UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes
ColumnPlot
ggsave(ColumnPlot,
       file = "Regulons_W0.MTb.MetaGeneSets.ALL.pdf",
       path = "GSEA_Figures",
       width = 18, height = 30, units = "in")
