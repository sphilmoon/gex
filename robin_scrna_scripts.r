library(SoupX)
library(ggplot2)
library(cowplot)
library(SingleR)
library(celldex)
library(Seurat)
library(SummarizedExperiment)
library(scran)
library(patchwork)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationHub)
library(dplyr)
library(enrichplot)
library(SingleCellExperiment)
library("reshape2")
library(viridis)



#REMOVE BACKGROUND IN NAIVE------------------------------------------------------------------------------------------------------
naive_toc = Read10X(data.dir ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_naive/outs/filtered_feature_bc_matrix/")
naive_tod = Read10X(data.dir ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_naive/outs/raw_feature_bc_matrix/")
naive_sc = SoupChannel(naive_tod, naive_toc)
naive_sc = SoupChannel(naive_tod, naive_toc, calcSoupProfile = FALSE)
naive_sc = estimateSoup(naive_sc)
naive_metadata <-read.csv(file ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_naive/outs/analysis/clustering/graphclust/clusters.csv")
naive_sc = setClusters(naive_sc, setNames(naive_metadata$Cluster, rownames(naive_metadata)))
# Estimate rho
naive_sc = autoEstCont(naive_sc)
# Clean the data
naive_out = adjustCounts(naive_sc)


#REMOVE BACKGROUND IN INFECTED----------------------------------------------------------------------------------------------------
infected_toc = Read10X(data.dir ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_infected/outs/filtered_feature_bc_matrix/")
infected_tod = Read10X(data.dir ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_infected/outs/raw_feature_bc_matrix/")
infected_sc = SoupChannel(infected_tod, infected_toc)
infected_sc = SoupChannel(infected_tod, infected_toc, calcSoupProfile = FALSE)
infected_sc = estimateSoup(infected_sc)


# Now loading the cluster from above
infected_metadata <-read.csv(file ="/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_infected/outs/analysis/clustering/graphclust/clusters.csv")
infected_sc = setClusters(infected_sc, setNames(infected_metadata$Cluster, rownames(infected_metadata)))
# Estimate rho
infected_sc = autoEstCont(infected_sc)
# Clean the data
infected_out = adjustCounts(infected_sc)


# Set up NAIVE object----------------------------------------------------------------------------------------
naive <- CreateSeuratObject(counts = naive_out, project = "naive", min.cells = 5)
naive$group <- "naive"
naive <- subset(naive, subset = nFeature_RNA > 500)
naive <- NormalizeData(naive, verbose = FALSE)
naive <- FindVariableFeatures(naive, selection.method = "vst", nfeatures = 2000)


# Set up INFECTED object
infected <- CreateSeuratObject(counts = infected_out, project = "infected", min.cells = 5)
infected$group <- "infected"
infected <- subset(infected, subset = nFeature_RNA > 500)
infected <- NormalizeData(infected, verbose = FALSE)
infected <- FindVariableFeatures(infected, selection.method = "vst", nfeatures = 2000)


#QUALITY CONTROL NAIVE-----------------------------------------------------------------------------------
naive[["percent.MT"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 
naive <- subset(naive, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.MT < 10)
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

#QUALITY CONTROL INFECTED--------------------------------------------------------------------------------
infected[["percent.MT"]] <- PercentageFeatureSet(infected, pattern = "^mt-")
VlnPlot(infected, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(infected, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(infected, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
infected <- subset(infected, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.MT < 10)
VlnPlot(infected, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)


 #Perform Integration-------------------------------------------------------------------------------------
pbmc.anchors <- FindIntegrationAnchors(object.list = list(naive, infected), dims = 1:20)
pbmc <- IntegrateData(anchorset = pbmc.anchors, dims = 1:20)
DefaultAssay(pbmc) <- "integrated"
pbmc


# Run the standard workflow for visualization and clustering-----------------------------------------------
pbmc <- ScaleData(pbmc, verbose = FALSE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = FALSE)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

# t-SNE and Clustering
pbmc <- RunUMAP(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Visualization------------------------------------------------------------------------------------------
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "group")
p2 <- DimPlot(pbmc, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


#UMAP------------------------------------------------------------------------------------------------
DimPlot(pbmc, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")
DimPlot(pbmc, reduction = "umap", label = TRUE)


#SAVE
saveRDS(pbmc, file = "/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_integrated.rds")


#Load in the data
pbmc <- readRDS(file = "/home/bmrc/robin/T_brucei_scRNA/20190704/brucei_integrated.rds")
pbmc

#--------------Single RAnnotation-----------------------------------------------------------------------
ref <- ImmGenData()
ref

whole_annotation <- SingleR(test = SummarizedExperiment(pbmc), ref = ref, assay.type.test=1, 
			                                labels = ref$label.main)

whole_annotation
head(pbmc@meta.data)
# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
pbmc$SingleR.pruned.calls <- whole_annotation$pruned.labels
pbmc$SingleR.calls <- whole_annotation$labels


#Run UMAP------------------------------------------------------------------------------------------
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.5)
head(pbmc@meta.data)
pbmc

#Annotation Diagnostics----------------------------------------------------------------------------
pred.grun <- SingleR(test=as.SingleCellExperiment(pbmc), ref=ref, labels=ref$label.main, de.method="wilcox")
table(pred.grun$labels)
plotScoreHeatmap(pred.grun)
plotDeltaDistribution(pred.grun, ncol = 3)

all.markers <- metadata(pred.grun)$de.genes
sceG$labels <- pred.grun$labels


# Next, switch the identity class of all cells to reflect replicate ID
Idents(pbmc) <- "SingleR.calls"

# How many cells are in each cluster
table(Idents(pbmc))

# How many cells are in each cluster
pbmc_counts <- table(Idents(pbmc))
head(pbmc_counts)

pbmc_long <- as.data.frame(pbmc_counts)
head(pbmc_long)
colnames(pbmc_long) <- c("cell_type", "Cell_count")
head(pbmc_long)
rownames(pbmc_long) <- pbmc_long$x

pbmc_long

bar_counts <- ggplot(pbmc_long,            # Create ggplot2 plot scaled to 1.00
		                   aes(x = cell_type, y = Cell_count,
				                         fill = cell_type)) +
  geom_bar(stat = "identity", color ="black") + 
    #scale_y_continuous(labels = scales::percent_format()) + 
    scale_x_discrete(name ="Cell type") +
      theme(legend.position="right") +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

bar_counts                               # Draw ggplot2 plot scaled to 1.00


# Can I create a Seurat object of just the Plasma cells and B cells?
B_Plasma_cells <- subset((pbmc), idents = c("B cells", "B cells, pro"))

# How many cells are in each cluster
table(Idents(B_Plasma_cells))

B_Plasma_cells
head(B_Plasma_cells@meta.data)
# Next, switch the identity class of all cells to reflect replicate ID
Idents(B_Plasma_cells) <- "seurat_clusters"


results <- SingleR(test = as.SingleCellExperiment(B_Plasma_cells), ref = ref, assay.type.test=1,
		                      labels = ref$label.fine)
results

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
B_Plasma_cells$SingleR.pruned.calls_B <- results$pruned.labels
B_Plasma_cells$SingleR.calls_B <- results$labels
head(B_Plasma_cells@meta.data)

#Run UMAP
B_Plasma_cells <- RunUMAP(B_Plasma_cells, dims = 1:10)
DimPlot(B_Plasma_cells, reduction = "umap", label = TRUE)
DimPlot(B_Plasma_cells, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_Plasma_cells, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_Plasma_cells, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)


# Next, switch the identity class of all cells to reflect replicate ID
Idents(B_Plasma_cells) <- "SingleR.calls_B"

# How many cells are in each cluster
table(Idents(B_Plasma_cells))

# Can I create a Seurat object of just the Plasma cells and B cells?
B_cells <- subset((B_Plasma_cells), idents = c("B cells (B.CD19CONTROL)","B cells (B.Fo)","B cells (B.FO)","B cells (B.FrE)",
					                                                      "B cells (B.FRE)","B cells (B.FrF)","B cells (B.GC)","B cells (B.MEM)",
											                                                     "B cells (B.MZ)","B cells (B.T1)","B cells (B.T2)","B cells (B.T3)",
											                                                     "B cells (B1a)","B cells (B1A)","B cells (B1b)","B cells (preB.FrC)",
																	                                                    "B cells (preB.FrD)","B cells (preB.FRD)","B cells (proB.CLP)","B cells (proB.FrA)",
																	                                                    "B cells (proB.FrBC)"))
# How many cells are in each cluster
B_cell_counts <- table(Idents(B_cells), B_cells$group)
head(B_cell_counts)

data_long <- as.data.frame(B_cell_counts)    
colnames(data_long) <- c("B_cell", "Group", "Relative_abundance")

rownames(data_long) <- data_long$x

data_long

ggp <- ggplot(data_long,            # Create ggplot2 plot scaled to 1.00
	                    aes(x = Group,
				                  y = Relative_abundance,
						                    fill = B_cell)) +
  geom_bar(position = "fill", stat = "identity", color ="black") + 
    scale_y_continuous(labels = scales::percent_format()) + 
      scale_x_discrete(name ="Group") +
        theme(legend.position="right")

ggp                                 # Draw ggplot2 plot scaled to 1.00



#Save counts in CSV
write.csv(B_cell_counts, file = "/home/bmrc/robin/T_brucei_scRNA/20190704/2021.03.26/B_cell_counts.txt")


#Run UMAP
B_cells <- RunUMAP(B_cells, dims = 1:10)
DimPlot(B_cells, reduction = "umap", label = TRUE)
DimPlot(B_cells, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_cells, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_cells, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

#DOT PLOT MARKER GENES------------------------------------------------------------------------------------
cd_genes <- c("Cd19", "Cd93", "Ly6d", "Ebf1","Ms4a1","Itgam","Zbtb32","Zbtb20","Cr2","Ighd","Fcer2a","Aicda","Mki67","Ighm","Sdc1","Ezh2","Bcl6", "Cd1d1", "Cd9","Pax5","Ptprc")
DotPlot(object = B_cells, features = cd_genes) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

#B Cells Only
FeaturePlot(B_cells, features = c("Cd79a","Cd19","Pax5", "Ebf1"), pt.size = 0.3)
VlnPlot(B_cells, features = c("Cd79a","Cd19","Pax5", "Ebf1"), pt.size = 0.3)


#Transitional B cells
VlnPlot(B_cells, features = c("Ly6d","Ebf1","Ms4a1", "Vpreb3"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Ly6d","Ebf1","Ms4a1", "Vpreb3"), pt.size = 0.3)

#B1 B cells / Atypical Memory B cells
VlnPlot(B_cells, features = c("Itgam","Zbtb32","Zbtb20"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Itgam","Zbtb32","Zbtb20"), pt.size = 0.3)

#Marginal Zone B cells
VlnPlot(B_cells, features = c("Ly6d","Ebf1","Cr2"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Ly6d","Ebf1","Cr2", "Ms4a1"), pt.size = 0.3)

#Follicular B cells
VlnPlot(B_cells, features = c("Ighd","Fcer2a", "Cd19"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Ighd","Fcer2a", "Cd19", "Cr2"), pt.size = 0.3)

#Germinal center B cells
VlnPlot(B_cells, features = c("Aicda", "Mki67"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Aicda", "Mki67"), pt.size = 0.3)

#Plasma cells
VlnPlot(B_cells, features = c("Sdc1","Xbp1","Ighm"), pt.size = 0.3)
FeaturePlot(B_cells, features = c("Sdc1","Xbp1","Ighm"), pt.size = 0.3)


# DIFFERENTIAL GENE EXPRESSION -------------------------------------------------------------------------------------
B_cell.markers <- FindAllMarkers(B_cells)
head(B_cell.markers)
write.csv(B_cell.markers, file= "/home/bmrc/robin/T_brucei_scRNA/May_2021/B_cell.markers.csv")

EnhancedVolcano(B_cell.markers,
		                lab = rownames(B_cell.markers),
				                x = 'avg_log2FC',
				                y = 'p_val',
						                pCutoff = 10e-4,
						                FCcutoff = 1.5,
								                xlim = c(-5.5, 5.5),
								                pointSize = 2.0,
										                labSize = 3.0,
										                title = 'Naive vs Infected',
												                subtitle = NULL,
												                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-4',
														                legendLabSize = 14,
														                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
																                colAlpha = 0.5,)

# Gene Ontology Analysis ------------------------------------------------------------------------------------------
head(B_cell.markers)
top100 <- B_cell.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
top100
top100pval <- subset(top100, rowSums(top100[5] < 0.05) > 0)
top100pval
write.csv(top100pval, file="/home/bmrc/robin/T_brucei_scRNA/May_2021/top100pval.csv")


# we want the log2 fold change 
original_gene_list <- top100pval$avg_log2FC
head(original_gene_list)

# name the vector
names(original_gene_list) <- top100pval$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gene_list

# Exctract significant results (padj < 0.05)
sig_genes = subset(top100pval, p_val_adj < 0.05)
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

#to get all the genes of all B cells into a universe
all.genes <- rownames(B_cells)
all.genes

#Keytypes
#keytypes(org.Mm.eg.db)

# ID Translator from SYMBOL to ENTREZID--------------------------------------------

genes = bitr(genes, fromType = "SYMBOL",
	                  toType = "ENTREZID",
			               OrgDb = org.Mm.eg.db)
head(genes)


#-----------Over representation GO Analysis---------------


ego_BP <- enrichGO(gene         = genes$SYMBOL,
		                      OrgDb         = org.Mm.eg.db,
				                         keyType       = "SYMBOL",
				                         ont           = "BP",
							                    pAdjustMethod = "BH",
							                    universe = all.genes,
									                       pvalueCutoff  = 0.01,
									                       qvalueCutoff  = 0.05)

ego_MF <- enrichGO(gene         = genes$SYMBOL,
		                      OrgDb         = org.Mm.eg.db,
				                         keyType       = "SYMBOL",
				                         ont           = "MF",
							                    pAdjustMethod = "BH",
							                    universe = all.genes,
									                       pvalueCutoff  = 0.01,
									                       qvalueCutoff  = 0.05)


ego_CC <- enrichGO(gene         = genes$SYMBOL,
		                      OrgDb         = org.Mm.eg.db,
				                         keyType       = "SYMBOL",
				                         ont           = "CC",
							                    pAdjustMethod = "BH",
							                    universe = all.genes,
									                       pvalueCutoff  = 0.01,
									                       qvalueCutoff  = 0.05)

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


#GOclusterplot <- compareCluster(geneCluster = gene_list, fun = "enrichGO", OrgDb = "org.Mm.eg.db")
#dotplot(GOclusterplot)


#------------------GO Classification---------------------------------------

ggoBP <- groupGO(gene     = genes$ENTREZID,
		                  OrgDb    = org.Mm.eg.db,
				                   keyType  ="ENTREZID",
				                   ont      = "BP",
						                    level    = 3,
						                    readable = TRUE)


ggoCC <- groupGO(gene     = genes$ENTREZID,
		                  OrgDb    = org.Mm.eg.db,
				                   keyType  ="ENTREZID",
				                   ont      = "CC",
						                    level    = 3,
						                    readable = TRUE)

ggoMF <- groupGO(gene     = genes$ENTREZID,
		                  OrgDb    = org.Mm.eg.db,
				                   keyType  ="ENTREZID",
				                   ont      = "MF",
						                    level    = 3,
						                    readable = TRUE)

barplot(ggoBP, drop=TRUE, showCategory=30)
barplot(ggoCC, drop=TRUE, showCategory=30)
barplot(ggoMF, drop=TRUE, showCategory=30)

#KEGG-----------------------------------------------------------------------------------
search_kegg_organism('mmu', by='kegg_code')
mmu <- search_kegg_organism('Mus musculus', by='scientific_name')
dim(mmu)
head(mmu)

kk <- enrichKEGG(gene         = genes$ENTREZID,
		                  organism     = 'mmu',
				                   pvalueCutoff = 0.1,
				                   qvalueCutoff = 0.1,
						                    pAdjustMethod = "BH",
						                    keyType = "GID")

head(kk)

#browseKEGG(kk, 'mmu05143')

barplot(kk, 
	        drop = TRUE, 
		        showCategory = 10, 
		        font.size = 8)

barplot(kk, drop=TRUE, showCategory=50, color='pvalue')

dotplot(kk, showCategory=50, color='pvalue')


kk1 <- pairwise_termsim(kk)
emapplot(kk1, color='pvalue')

heatplot(kk, showCategory = 30, foldChange = )

cnetplot(kk, categorySize="pvalue", foldChange = gene_list, circular=TRUE)

write.csv(kk, file="/home/bmrc/robin/T_brucei_scRNA/May_2021/KEGG2 UP enriched.csv")


kk2 <- gseKEGG(geneList     = genes$ENTREZID,
	                      organism     = 'mmu',
			                     nPerm        = 1000,
			                     minGSSize    = 120,
					                    pvalueCutoff = 0.05,
					                    verbose      = FALSE)
head(kk2)





