# GSEA with CluterProfiler
# E. Lamont
# 4/21/25

### *** Won't work because they don't have Mtb as the organism.....

# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

################################################
################ LOAD PACKAGES #################

BiocManager::install("clusterProfiler", version = "3.20")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)

################################################
################# LOAD DATA ####################

# This requires log2fold change data... from DEG...

# Use W2.ComparedTo.W0

################################################
############### PREPARE INPUT ##################

# reading in data from deseq2
# df = read.csv("drosphila_example_de.csv", header=TRUE)
df <- W2.ComparedTo.W0

# we want the log2 fold change 
# original_gene_list <- df$log2FoldChange
original_gene_list <- df$LOG2FOLD

# name the vector
# names(original_gene_list) <- df$X
names(original_gene_list) <- df$GENE_ID

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


################################################
########### GENE SET ENRICHMENT ################

gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

### Won't work because they don't have Mtb as the organism.....





