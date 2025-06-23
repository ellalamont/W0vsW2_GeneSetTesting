# Looking at the MetaGeneSets from Bob's pipeline and subsetting after everything has been run - W0 vs broth!!
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
############## iModulons: JUST CENTRAL CARBON #############

CentralCarbon_iModulons <- c("Peptidoglycan Biosynthesis", "Central Carbon Metabolism", "Fumarate Reductase", "PrpR", "BkaR", "Nicotinate Metabolism") # Needs the \\ to detect it
pattern <- str_c(CentralCarbon_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_CentralCarbon <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_CentralCarbon %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  # arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.78, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.8, 3.2), expand = c(0, 0)) +
  labs(x="Central Carbon iModulon", y="LOG2FOLD",
       title="Central Carbon iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_CentralCarbon.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
######### iModulons: JUST AMINO ACID BIOSYNTHESIS #########

AminoAcid_iModulons <- c("GroEL-GroES Complex", "Leucine Related", "LysG", "ArgR") # Needs the \\ to detect it
pattern <- str_c(AminoAcid_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_AminoAcid <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_AminoAcid %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
  labs(x="Amino Acid Biosynthesis iModulon", y="LOG2FOLD",
       title="Amino Acid Biosynthesis iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_AminoAcid.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
############## iModulons: JUST NUCLEIC ACIDS ##############

NucleicAcid_iModulons <- c("PyrR", "Rv0135\\+Rv1019", "Nucleic Acid Hydrolysis") # Needs the \\ to detect it
pattern <- str_c(NucleicAcid_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_NucleiAcid <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_NucleiAcid %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-2.5, 1), expand = c(0, 0)) +
  labs(x="Nucleic Acid iModulon", y="LOG2FOLD",
       title="Nucleic Acid iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_NucleicAcid.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
######## iModulons: JUST FATTY ACID/CHOLESTEROL ###########

FattyAcid.Cholesterol_iModulons <- c("Fatty Acid Biosynthesis", "KstR2", "Mycofactocin Synthesis Pathway", "FasR", "Polyketide Synthase Complex", "Rv0681") # Needs the \\ to detect it
pattern <- str_c(FattyAcid.Cholesterol_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_FattyAcid.Cholesterol <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_FattyAcid.Cholesterol %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-2.5, 3.5), expand = c(0, 0)) +
  labs(x="FattyAcid/Cholesterol iModulon", y="LOG2FOLD",
       title="FattyAcid/Cholesterol iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_FattyAcid.Cholesterol.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
############# iModulons: JUST METAL RELATED ###############

Metal_iModulons <- c("RicR", "IdeR", "M-box", "Zur", "Hpt-2b Induced") # Needs the \\ to detect it
pattern <- str_c(Metal_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_Metal <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Metal %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0, 0)) +
  labs(x="Metal iModulon", y="LOG2FOLD",
       title="Metal iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Metal.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
########### iModulons: JUST SULFUR METABOLISM #############

SulfurMetabolism_iModulons <- c("Sulfur Metabolism") # Needs the \\ to detect it
pattern <- str_c(SulfurMetabolism_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_SulfurMetabolism <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_SulfurMetabolism %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.5, 2.5), expand = c(0, 0)) +
  labs(x="Sulfur Metabolism iModulon", y="LOG2FOLD",
       title="Sulfur Metabolism iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_SulfurMetabolism.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 3, units = "in")


###########################################################
################# iModulons: JUST GROWTH ##################

Growth_iModulons <- c("Positive Regulation of Growth") # Needs the \\ to detect it
pattern <- str_c(Growth_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_Growth <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Growth %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
  labs(x="Growth iModulon", y="LOG2FOLD",
       title="Growth iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Growth.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 3, units = "in")


###########################################################
################ iModulons: JUST REDOX ####################

Redox_iModulons <- c("DevR-1", "WhiB4", "DevR-2", "WhiB1", "WhiB4/IdeR", "Rv1828/SigH", "Rv1776c\\+WhiB4", "VirS", "WhiB6") # Needs the \\ to detect it
pattern <- str_c(Redox_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_Redox <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Redox %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -2.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-2.5, 0.5), expand = c(0, 0)) +
  labs(x="Redox iModulon", y="LOG2FOLD",
       title="Redox iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Redox.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")

###########################################################
############## iModulons: JUST ACID STRESS ################

AcidStress_iModulons <- c("MarR") # Needs the \\ to detect it
pattern <- str_c(AcidStress_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_AcidStress <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_AcidStress %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.48, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-0.5, 3), expand = c(0, 0)) +
  labs(x="Acid Stress iModulon", y="LOG2FOLD",
       title="Acid Stress iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_AcidStress.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 3, units = "in")


###########################################################
############## iModulons: JUST ANTIBIOTIC #################

# Lsr2 not present in my data!!
Antibiotic_iModulons <- c("Lsr2", "Blal", "Rv0078\\+Rv2034", "WhiB7", "IniR") # Needs the \\ to detect it
pattern <- str_c(Antibiotic_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_Antibiotic <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Antibiotic %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -2.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-3, 3), expand = c(0, 0)) +
  labs(x="Antibiotic iModulon", y="LOG2FOLD",
       title="Antibiotic iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Antibiotic.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")


###########################################################
######### iModulons: JUST VIRULENCE/PERSISTENCE ###########

# Lsr2 not present in my data!!
Virulence.Persistence_iModulons <- c("Rv0576", "Mce1R", "SigH", "PhoP", "Mce3R", "MprA", "PDIM\\;PGL Synthesis", "Rv2488c", "SigC", "SigD", "MbcA\\+Rv3249c\\+Rv3066", "SigK") # Needs the \\ to detect it
pattern <- str_c(Virulence.Persistence_iModulons, collapse = "|") # Collapse all the things I am interested in into a pattern separated by or

# Filter
MetaGeneSets_W0vsBroth_ALL_iModulons_Virulence.Persistence <- MetaGeneSets_W0vsBroth_ALL_iModulons %>% filter(str_detect(PathName, pattern))

ColumnPlot <- MetaGeneSets_W0vsBroth_ALL_iModulons_Virulence.Persistence %>% 
  mutate(PathName = str_wrap(PathName, width = 60)) %>%  # Adjust width as needed
  arrange(desc(LOG2FOLD)) %>%
  ggplot(aes(reorder(PathName, LOG2FOLD), LOG2FOLD)) + 
  geom_col(aes(fill = AVG_PVALUE < 0.05)) + 
  scale_fill_manual(values=c("FALSE" = "#999999", "TRUE" = "red3")) + 
  coord_flip() +
  geom_text(aes(y = -0.98, label = paste0("n = ", N_Genes)), hjust = 0, size = 2.5) +
  scale_y_continuous(limits = c(-1, 3.5), expand = c(0, 0)) +
  labs(x="Virulence/Persistence iModulon", y="LOG2FOLD",
       title="Virulence/Persistence iModulons in W0 Sputum vs Log broth \nusing UP and DOWN combined",
       subtitle = "The UP with LOG2FOLD > 0 and the DOWN with LOG2FOLD < 0") + 
  my_plot_themes + theme(axis.text.y = element_text(size=12))
ColumnPlot
ggsave(ColumnPlot,
       file = "iModulons_W0.MTb.MetaGeneSets.ALL_Virulence.Persistence.pdf",
       path = "GSEA_Figures/MetaGeneSets/iModulons/W0vsBroth",
       width = 10, height = 6, units = "in")




