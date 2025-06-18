# Looking at the MetaGeneSets from Bob's pipeline and subsetting after everything has been run
# E. Lamont
# 6/18/25

# Following from MetaGeneSets_Analysis.R, but just looking at iModulons and their categories

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
############## iModulons: JUST THE VIRULENCE ##############
# Just plot all the virulence/persistence associated iModulons

# Rv0576, Mce1R, SigH, PhoP, Mce3R, MprA, PDIM;PGL Synthesis, Rv2488c, SigC, SigD, MbcA+Rv3249c+Rv3066
Virulence.Persistence_iModulons <- c("Rv0576", "Mce1R", "SigH", "PhoP", "Mce3R", "MprA", "PDIM;PGL Synthesis", "Rv2488c", "SigC", "SigD", "MbcA\\+Rv3249c\\+Rv3066") # Needs the \\ to detect it
pattern <- str_c(Virulence.Persistence_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter to have just the virulence/persistence iModulons
MetaGeneSets_W0vsBroth_ALL_iModulons_Virulence.Persistence <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Virulence.Persistence %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.78, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.8, 3.5), expand = c(0, 0)) +
  labs(x="Virulence/Persistence iModulon", y="LOG2FOLD",
       title="Virulence/Persistence iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Virulence.Persistence.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons",
       width = 10, height = 6, units = "in")





