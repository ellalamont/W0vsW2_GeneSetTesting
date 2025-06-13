# Trying fgsea
# E. Lamont 
# 6/12/15

# https://github.com/alserglab/fgsea
# https://bioconductor.org/packages/devel/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
# https://biostatsquid.com/fgsea-tutorial-gsea/
# https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package

source("Import_data.R")

library(fgsea)
library(data.table)
library(ggplot2)

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

set.seed(42) # The p-values have some random changes going here...
GSEAres <- fgsea(pathways = my_geneSet, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 2,
                 maxSize = 500,
                 nproc = 1) # for parallelisation


# Top 6 enriched pathways (ordered by p-val)
head(GSEAres[order(pval), ])

sum(GSEAres[, padj < 0.05]) # 9 significant pathways only when adjusting p-values
sum(GSEAres[, pval < 0.05]) # 14 if not adjusted

topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = 4), pathway]
topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = 4), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)
plotGseaTable(my_geneSet[topPathways], stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()


# https://stephenturner.github.io/deseq-to-fgsea/#using_the_fgsea_package

set.seed(42)
fgseaRes <- fgsea(pathways=my_geneSet, stats=rankings)
# Bascially the same as above (if set.seed same would be the same)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()


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
       file = "iModulons_fgsea_ColPlot_4.pdf",
       path = "GSEA_Figures",
       width = 18, height = 18, units = "in")


# Just visualize the ones with NES > |1|
GSEA_barPlot <- GSEAres %>% 
  filter(abs(NES) >= 1) %>%
  arrange(desc(NES)) %>% 
  ggplot(aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  scale_fill_manual(values=c("#999999", "red3")) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", size)), hjust = 0, size = 2.5) +
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="Normalized Enrichment Score",
       title="iModulons in W0 Sputum vs Log broth using fgsea",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "iModulons_fgsea_ColPlot_subset.pdf",
       path = "GSEA_Figures",
       width = 18, height = 10, units = "in")


###########################################################
############### Mtb.Regulons: RUN THE FGSEA ###############

set.seed(42) # The p-values have some random changes going here...
GSEAres_Regulons <- fgsea(pathways = allGeneSetList$MTb.Regulons, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 2,
                 maxSize = 500,
                 nproc = 1) # for parallelisation

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
       file = "Regulons_fgsea_ColPlot_4.pdf",
       path = "GSEA_Figures",
       width = 18, height = 30, units = "in")






