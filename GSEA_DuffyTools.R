# Trying Bob's GSEA functions
# E. Lamont 
# 5/2/25

# TO DO:
# Pipe.GeneSetAnalysis from Duffy Tools: https://github.com/robertdouglasmorrison/DuffyTools/blob/master/man/pipe.GeneSetAnalysis.Rd
# Look at what is currently run on the lenovo and see if I can give it specific gene sets


# Stop scientific notation
options(scipen = 999) 
# options(scipen = 0) # To revert back to default


# https://rdrr.io/github/robertdouglasmorrison/DuffyTools/src/R/pipe.GSEA.R


do.GSEA(group1 = "W0", prefix = "MTb", path="/Users/snork-maiden/Documents/Micro_grad_school/Sherman_Lab/R_projects/Sputum/W0vsW2_GeneSetTesting/JOINED_BobAverages/MTb.MetaResults.W0_vs_Broth", tool = "MetaResults")


# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

BiocManager::install("fgsea")
library(fgsea)
library(data.table)

data(examplePathways)
data(exampleRanks)
set.seed(42)


# https://biostatsquid.com/gene-set-enrichment-analysis/
# https://biostatsquid.com/fgsea-tutorial-gsea/
# For GSEA, genes need to be ranks

# Duffy Tools ranking function
# diffExpressRankOrder
# rank_test <- diffExpressRankOrder(folds = W0.ComparedTo.Broth$LOG2FOLD, pvalues = W0.ComparedTo.Broth$AVG_PVALUE)
# This is ranking starting at 1, nothing negative which I think is what I need


###########################################################
################# FGSEA SQUIDTIPS METHOD ##################

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

## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = allGeneSetList$Ella_GeneSets, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 5,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
# The p-values change every time I run this......

# Not working
# topPathwaysUp <- GSEAres[ES > 0][head(order(pval), n = number_of_top_pathways_up), pathway]
# topPathwaysDown <- GSEAres[ES < 0][head(order(pval), n = number_of_top_pathways_down), pathway]
# topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
#pdf(file = paste0(filename, '_gsea_top30pathways.pdf'), width = 20, height = 15)

plotGseaTable(allGeneSetList$Ella_GeneSets, stats = rankings, fgseaRes = GSEAres, gseaParam = 0.5)
#dev.off()


###########################################################
################# FGSEA POSITIVE RANKINGS #################
# Use the Duffy tools ranking system then switch to positive when running the fgsea to see what happens

Rankings_Duffy <- diffExpressRankOrder(folds = W0.ComparedTo.Broth$LOG2FOLD, pvalues = W0.ComparedTo.Broth$AVG_PVALUE)
names(Rankings_Duffy) <- W0.ComparedTo.Broth$GENE_ID # genes as names#

GSEAres <- fgsea(pathways = allGeneSetList$Ella_GeneSets, # List of gene sets to check
                 stats = Rankings_Duffy,
                 scoreType = 'pos', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 minSize = 5,
                 maxSize = 500,
                 nproc = 1) # for parallelisation
# The p-values are really bad now! 



