# Gene Ontology analysis (CD4)

# load the library packages
library(monocle3)
library(ggplot2)
library(cowplot)
library(Seurat)
library(scran)
library(patchwork)
library(clusterProfiler)
library(dplyr)
library(enrichplot)
library(org.Mm.eg.db)
library(AnnotationHub)
library(SingleCellExperiment)
library("reshape2")
library(viridis)

# open the cd4 gex file
cd4_gex_073121 <- read.csv(file= "/home/bmrc/Desktop/cd4_gex_073121.csv")
head(cd4_gex_073121)

# find top 100 genes
cd4_top100 <- cd4_gex_073121 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# sort based on p < 0.05
cd4_top100_pval <- subset(cd4_top100, rowSums(cd4_top100[5] < 0.05) > 0)

# we want the log2 fold change 
original_gene_list <- cd4_top100_pval$avg_log2FC
head(original_gene_list)

# name the vector
names(original_gene_list) <- cd4_top100_pval$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list

# Exctract significant results (padj < 0.05)
sig_genes = subset(cd4_top100_pval, p_val_adj < 0.05)
sig_genes

# From significant results, we want to filter on log2fold change
genes <- sig_genes$avg_log2FC

# Name the vector
names(genes) <- sig_genes$gene

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]
genes

# to get all the genes of all B cells into a universe
all.genes <- rownames(B_cells)
all.genes

# Keytypes
# keytypes(org.Mm.eg.db)

# ID Translator from SYMBOL to ENTREZID--------------------------------------------

genes = bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
head(genes)


#-----------Over representation GO Analysis---------------


ego_BP <- enrichGO(gene = genes$SYMBOL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   universe = all.genes,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

ego_MF <- enrichGO(gene = genes$SYMBOL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   universe = all.genes,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)


ego_CC <- enrichGO(gene = genes$SYMBOL,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   universe = all.genes,
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)

barplot(ego_BP,showCategory=20, color='pvalue')
barplot(ego_MF,showCategory=20, color='pvalue')
barplot(ego_CC,showCategory=20, color='pvalue')

dotplot(ego_BP)
dotplot(ego_MF)
dotplot(ego_CC)

cnetplot(ego_BP, categorySize="pvalue", foldChange = gene_list, circular=FALSE)
cnetplot(ego_BP, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_MF, categorySize="pvalue", foldChange = gene_list, circular=FALSE)
cnetplot(ego_MF, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_CC, categorySize="pvalue", foldChange = gene_list, circular=FALSE)
cnetplot(ego_CC, foldChange=gene_list, circular = TRUE, colorEdge = TRUE)

heatplot(ego_BP, showCategory = 20, foldChange = gene_list)
heatplot(ego_MF, showCategory = 20, foldChange = gene_list)
heatplot(ego_CC, showCategory = 20, foldChange = gene_list)

ego_BP1 <- pairwise_termsim(ego_BP)
emapplot(ego_BP1, color='pvalue')




#------------------GO Classification---------------------------------------

ggoBP <- groupGO(gene = genes$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 keyType ="ENTREZID",
                 ont = "BP",
                 level = 3,
                 readable = TRUE)


ggoCC <- groupGO(gene = genes$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 keyType ="ENTREZID",
                 ont = "CC",
                 level = 3,
                 readable = TRUE)

ggoMF <- groupGO(gene = genes$ENTREZID,
                 OrgDb = org.Mm.eg.db,
                 keyType ="ENTREZID",
                 ont = "MF",
                 level = 3,
                 readable = TRUE)

barplot(ggoBP, drop=TRUE, showCategory=30)
barplot(ggoCC, drop=TRUE, showCategory=30)
barplot(ggoMF, drop=TRUE, showCategory=30)