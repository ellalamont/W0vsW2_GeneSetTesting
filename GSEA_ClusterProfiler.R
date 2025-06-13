# GSEA with CluterProfiler
# E. Lamont
# 6/13/25

# Got ChatGPT to help me with this one

# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)

# Prepare gene ranks
rankings <- sign(W0.ComparedTo.Broth$LOG2FOLD)*(-log10(W0.ComparedTo.Broth$AVG_PVALUE))
names(rankings) <- W0.ComparedTo.Broth$GENE_ID # genes as names#
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
geneList <- rankings

# Convert list to data.frame for clusterProfiler
geneset_df <- stack(allGeneSetList$MTb.iModulons)
colnames(geneset_df) <- c("gene", "set")  # must be columns named "gene" and "set"

geneset_df <- stack(allGeneSetList$MTb.iModulons)
colnames(geneset_df) <- c("gene", "term")  # term = pathway name
geneset_df <- geneset_df[, c("term", "gene")]  # TERM2GENE = (term, gene)

set.seed(42)
gsea_res <- GSEA(
  geneList = geneList,
  TERM2GENE = geneset_df,
  verbose = FALSE,
  minGSSize = 2,
  maxGSSize = 500,
  pvalueCutoff = 1
)
head(as.data.frame(gsea_res))

# Top plot
ridgeplot(gsea_res, showCategory = 10)

fgseaResTidy2 <- gsea_res %>%
  as_tibble() %>%
  arrange(desc(NES))

ggplot(fgseaResTidy2, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill=p.adjust<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

upsetplot(gsea_res) 



top_terms <- gsea_res@result %>%
  as_tibble() %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = 10) %>%
  mutate(term = reorder(Description, NES))

ggplot(top_terms, aes(x = NES, y = term, fill = p.adjust)) +
  geom_col() +
  scale_fill_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    title = "Top Enriched Pathways"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))


###########################################################
################## iModulon COLUMN PLOT ###################
# Haven't gotten this to work yet



GSEA_barPlot <- gsea_res %>% 
  arrange(desc(NES)) %>%
  ggplot(aes(reorder(Description, NES), NES)) + 
  geom_col(aes(fill = p.adjust<0.05)) + 
  geom_text(aes(y = -1.98, label = paste0("n = ", setSize)), hjust = 0, size = 2.5) +
  scale_fill_manual(values=c("#999999", "red3")) + 
  coord_flip() +
  scale_y_continuous(limits = c(-2, 2), expand = c(0, 0)) +
  labs(x="iModulon", y="Normalized Enrichment Score",
       title="iModulons in W0 Sputum vs Log broth using ClusterProfiler",
       subtitle = "min gene set size = 2") + 
  my_plot_themes
GSEA_barPlot
ggsave(GSEA_barPlot,
       file = "iModulons_fgsea_ColPlot_4.pdf",
       path = "GSEA_Figures",
       width = 18, height = 18, units = "in")




