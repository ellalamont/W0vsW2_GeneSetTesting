# Trying fgsea
# E. Lamont 
# 6/12/15

# https://github.com/alserglab/fgsea
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://biostatsquid.com/fgsea-tutorial-gsea/
# https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package


source("Import_data.R")

library(fgsea)
# library(data.table)

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
##################### RANK THE GENES ######################

# https://biostatsquid.com/fgsea-tutorial-gsea/
# Using the squidtip method
rankings <- sign(W0.ComparedTo.Broth$LOG2FOLD)*(-log10(W0.ComparedTo.Broth$AVG_PVALUE))
names(rankings) <- W0.ComparedTo.Broth$GENE_ID # genes as names#
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)

# Check that there aren't any inf numbers
max(rankings)
min(rankings)

# Check the rankings of genes
ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###########################################################
################# iMODULONS: RUN THE FGSEA ################

# allGeneSetList$MTb.iModulons

my_geneSet <- allGeneSetList$MTb.iModulons

set.seed(23) # The p-values have some random changes going here...
GSEAres <- fgsea(pathways = my_geneSet, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 2,
                 maxSize = 500,
                 nproc = 1) # for parallelisation


# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

sum(GSEAres[, padj < 0.05]) # 9 significant pathways only when adjusting p-values
sum(GSEAres[, pval < 0.05]) # 18 if not adjusted

topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = 4), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = 4), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(my_geneSet[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()



###########################################################
################## iModulon COLUMN PLOT ###################

# Using GSEAres although they are the same

GSEA_barPlot <- GSEAres %>% 
  arrange(desc(NES)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="Normalized Enrichment Score",
       title="iModulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "iModulons_W0vsBroth_fgsea.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 18, units = "in")


# Just the significant ones
GSEA_barPlot <- GSEAres %>% 
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  scale_fill_manual(values=c("red3")) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="Normalized Enrichment Score",
       title="iModulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "iModulons_W0vsBroth_fgsea_SignificantOnly.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 7, units = "in")


###########################################################
############### Mtb.Regulons: FGSEA + PLOT ################

set.seed(23) # The p-values have some random changes going here...
GSEAres_Regulons <- fgsea(pathways = allGeneSetList$MTb.Regulons, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 2,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

sum(GSEAres_Regulons[, padj < 0.05]) # 7 significant pathways only when adjusting p-values
sum(GSEAres_Regulons[, pval < 0.05]) # 30 if not adjusted

GSEA_barPlot <- GSEAres_Regulons %>% 
  arrange(desc(NES)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Regulons", y="Normalized Enrichment Score",
       title="Regulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "Regulons_W0vsBroth_fgsea.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 30, units = "in")

# Just the significant ones
GSEA_barPlot <- GSEAres_Regulons %>% 
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  scale_fill_manual(values=c("red3")) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Regulon", y="Normalized Enrichment Score",
       title="Regulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "Regulons_W0vsBroth_fgsea_SignificantOnly.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 7, units = "in")


###########################################################
############ Mtb.TFOE.Regulons: FGSEA + PLOT ##############

set.seed(42) # The p-values have some random changes going here...
GSEAres_TFOE.Regulons <- fgsea(pathways = allGeneSetList$MTb.TFOE.Regulons, # List of gene sets to check
                          stats = rankings,
                          scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                          minSize = 2,
                          maxSize = 500,
                          nproc = 1) # for parallelisation

sum(GSEAres_TFOE.Regulons[, padj < 0.05]) # 11 significant pathways only when adjusting p-values
sum(GSEAres_TFOE.Regulons[, pval < 0.05]) # 31 if not adjusted

GSEA_barPlot <- GSEAres_TFOE.Regulons %>% 
  arrange(desc(NES)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="TFOE.Regulons", y="Normalized Enrichment Score",
       title="TFOE.Regulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "TFOE.Regulons_W0vsBroth_fgsea.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 30, units = "in")


# Just the significant ones
GSEA_barPlot <- GSEAres_TFOE.Regulons %>% 
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  scale_fill_manual(values=c("red3")) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="TFOE.Regulon", y="Normalized Enrichment Score",
       title="TFOE.Regulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "TFOE.Regulons_W0vsBroth_fgsea_SignificantOnly.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 7, units = "in")


###########################################################
############ Tuberculist.GO.Ontology: FGSEA + PLOT ##############

set.seed(23) # The p-values have some random changes going here...
GSEAres_Tuberculist.GO.Ontology <- fgsea(pathways = allGeneSetList$MTb.Tuberculist.GO.Ontology, # List of gene sets to check
                               stats = rankings,
                               scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                               minSize = 2,
                               maxSize = 500,
                               nproc = 1) # for parallelisation

sum(GSEAres_Tuberculist.GO.Ontology[, padj < 0.05]) # 8 significant pathways only when adjusting p-values
sum(GSEAres_Tuberculist.GO.Ontology[, pval < 0.05]) # 47 if not adjusted

GSEA_barPlot <- GSEAres_Tuberculist.GO.Ontology %>% 
  arrange(desc(NES)) %>%
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.Ontology", y="Normalized Enrichment Score",
       title="Tuberculist.GO.Ontology in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "Tuberculist.GO.Ontology_W0vsBroth_fgsea.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 49.9, units = "in")


# Just the significant ones
GSEA_barPlot <- GSEAres_Tuberculist.GO.Ontology %>% 
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  scale_fill_manual(values=c("red3")) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="Tuberculist.GO.Ontology", y="Normalized Enrichment Score",
       title="Tuberculist.GO.Ontology in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "Tuberculist.GO.Ontology_W0vsBroth_fgsea_SignificantOnly.pdf",
       path = "GSEA_Figures/fgsea_package",
       width = 18, height = 7, units = "in")


