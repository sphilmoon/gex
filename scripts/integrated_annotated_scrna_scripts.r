library(dplyr)
library(Seurat)
library(patchwork)
library(SoupX)
library(ggplot2)
library(celldex)
library(SingleR)
library(org.Mm.eg.db)
library(SingleCellExperiment)
library(SummarizedExperiment)
# library(pheatmap)

#########################################################################################################
## naive
#########################################################################################################

# 1. SoupX (pre-processing) -----------------------------------------------------------------------------

# Load the the naive data
naive_toc <- Read10X(data.dir = "/home/phil/naive/outs/filtered_feature_bc_matrix/")
naive_tod <- Read10X(data.dir = "/home/phil/naive/outs/raw_feature_bc_matrix/")
naive_sc <- SoupChannel(naive_tod, naive_toc, calcSoupProfile = FALSE)
naive_sc <- estimateSoup(naive_sc)

# Now loading the cluster from above
naive_cluster <-read.csv(file ="/home/phil/naive/outs/analysis/clustering/graphclust/clusters.csv")
naive_sc = setClusters(naive_sc, setNames(naive_cluster$Cluster, rownames(naive_cluster)))

# estimating the contamination fraction
naive_sc = autoEstCont(naive_sc)
# removing the contamination, then yield a new soupX variable
naive_soupx = adjustCounts(naive_sc)


# 2. Clustering (Seurat) -----------------------------------------------------------------------------

# a UMI count matrix output of the 'cellranger count' pipeline from 10X.
# contains number of molecules for each feature (row) in each cell (column).


# Initialize the Seurat object with the raw (non-normalized data).
naive_data <- CreateSeuratObject(counts = naive_soupx, project = "naive", min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
naive_data[["percent.mt"]] <- PercentageFeatureSet(naive_data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(naive_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(naive_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(naive_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Plot 2 is working in this case R^2 0.87 

# Normalization 
naive_data <- subset(naive_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
naive_data <- NormalizeData(naive_data, normalization.method = "LogNormalize", scale.factor = 10000)
naive_data <- NormalizeData(naive_data)

naive_data <- FindVariableFeatures(naive_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(naive_data), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(naive_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#########################################################################################################
## 14 dpi
#########################################################################################################

# 1. SoupX (pre-processing) -----------------------------------------------------------------------------

# Load the the 14 dpi data
infected_toc <- Read10X(data.dir = "/home/phil/infected/outs/filtered_feature_bc_matrix/")
infected_tod <- Read10X(data.dir = "/home/phil/infected/outs/raw_feature_bc_matrix/")
infected_sc <- SoupChannel(infected_tod, infected_toc, calcSoupProfile = FALSE)
infected_sc <- estimateSoup(infected_sc)

# Now loading the cluster from above
infected_cluster <-read.csv(file ="/home/phil/infected/outs/analysis/clustering/graphclust/clusters.csv")
infected_sc = setClusters(infected_sc, setNames(infected_cluster$Cluster, rownames(infected_cluster)))

# estimating the contamination fraction
infected_sc = autoEstCont(infected_sc)
# removing the contamination, then yield an output file
infected_soupx = adjustCounts(infected_sc)

# 2. Clustering (Seurat) -----------------------------------------------------------------------------

# a UMI count matrix output of the 'cellranger count' pipeline from 10X.
# contains number of molecules for each feature (row) in each cell (column).


# Initialize the Seurat object with the raw (non-normalized data).
infected_data <- CreateSeuratObject(counts = infected_soupx, project = "infected", min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
infected_data[["percent.mt"]] <- PercentageFeatureSet(infected_data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(infected_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(infected_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(infected_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
# Plot 2 is working in this case R^2 0.87 

# Normalization 
infected_data <- subset(infected_data, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)
infected_data <- NormalizeData(infected_data, normalization.method = "LogNormalize", scale.factor = 10000)
infected_data <- NormalizeData(infected_data)

infected_data <- FindVariableFeatures(infected_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(infected_data), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(infected_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


#########################################################################################################
## Integrating naive and 14 dpi datasets 
#########################################################################################################

# Perform Integration-------------------------------------------------------------------------------------
integrated_data.anchors <- FindIntegrationAnchors(object.list = list(naive_data, infected_data), dims = 1:20)
integrated_data <- IntegrateData(anchorset = integrated_data.anchors, dims = 1:20)
DefaultAssay(integrated_data) <- "integrated"
integrated_data

all.genes <- rownames(integrated_data)
integrated_data <- ScaleData(integrated_data, features = all.genes)
integrated_data <- RunPCA(integrated_data, features = VariableFeatures(object = integrated_data))
# Examine and visualize PCA results a few different ways
print(integrated_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(integrated_data, dims = 1:2, reduction = "pca")
DimPlot(integrated_data, reduction = "pca")

# PC1 - PC15
DimHeatmap(integrated_data, dims = 1:15, cells = 500, balanced = TRUE)

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
integrated_data <- JackStraw(integrated_data, num.replicate = 100)
integrated_data <- ScoreJackStraw(integrated_data, dims = 1:20)
JackStrawPlot(integrated_data, dims = 1:15)
ElbowPlot(integrated_data)

# tSNE and clustering
integrated_data <- FindNeighbors(integrated_data, dims = 1:20)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(integrated_data), 5)

### Running non-linear dimenstioan reduction (UMAP or tSNE) ###

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
integrated_data <- RunUMAP(integrated_data, dims = 1:20)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(integrated_data, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(integrated_data, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated_data, reduction = "umap", label = TRUE)

### I can save umap, the non-linear dimensional reduction at this point. ###
saveRDS(integrated_data, file = "clustering_integrated.rds")


"
#########################################################################################################
## Robin's integration 
#########################################################################################################

# 3. Integration-------------------------------------------------------------------------------------
integrated.anchors <- FindIntegrationAnchors(object.list = list(naive_data, infected_data), dims = 1:20)
integrated <- IntegrateData(anchorset = integrated.anchors, dims = 1:20)
DefaultAssay(integrated) <- "integrated"
integrated


# Run the standard workflow for visualization and clustering-----------------------------------------------
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
VizDimLoadings(integrated, dims = 1:2, reduction = "pca")
DimPlot(integrated, reduction = "pca")
DimHeatmap(integrated, dims = 1, cells = 500, balanced = TRUE)

# t-SNE and Clustering
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:20)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.5)

# Visualization------------------------------------------------------------------------------------------
p1 <- DimPlot(integrated, reduction = "umap", group.by = "group")
p2 <- DimPlot(integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


#UMAP------------------------------------------------------------------------------------------------
DimPlot(integrated, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated, reduction = "umap", label = TRUE)

#SAVE
saveRDS(integrated, file = "/home/phil/integrated_brucei.rds")
"


"
### // OPTIONAL //
### Finding differentially expressed features (e.g. cluster biomarkers) ###

# find all markers of cluster 2
cluster2.markers <- FindMarkers(infected_data, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(infected_data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
# this process will take some time
infected_data.markers <- FindAllMarkers(infected_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
infected_data.markers %>%
	    group_by(cluster) %>%
	        top_n(n = 2, wt = avg_log2FC)

	cluster0.markers <- FindMarkers(infected_data, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
	VlnPlot(infected_data, features = c("Ngp", "Gm26917"))

	# plottin raw counts as well
	VlnPlot(infected_data, features = c("Ngp", "Gm26917"), slot = "counts", log = TRUE)

	# plotting multiple feature plots
	FeaturePlot(infected_data, features = c("Ngp", "Gm26917", "Ppib", "C1qa", "Cd3g", "Apol11b"))

	# plot a heatamap
	infected_data.markers %>%
		    group_by(cluster) %>%
		        top_n(n = 10, wt = avg_log2FC) -> top10
		DoHeatmap(infected_data, features = top10$gene) + NoLegend()

		# save the rds file before annotation step
		saveRDS(infected_data, file = "infected_umap_trial.rds")
		"

		# 3. Annotation using SingleR package -----------------------------------------------------------------------------

		integrated <- readRDS(file = "/home/phil/integrated_brucei/clustering_integrated.rds")
		integrated

		# load ImmGenData library
		ref <- ImmGenData()
		ref

		"
		# annotate integrated variable
		whole_annotation <- SingleR(test = integrated, ref = ref, assay.type.test=1, labels = ref$label.main, drop = FALSE)
		whole_annotation
		head(integrated@meta.data)

		# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
		integrated$SingleR.pruned.calls <- whole_annotation$pruned.labels
		integrated$SingleR.calls <- whole_annotation$labels

		#Run UMAP------------------------------------------------------------------------------------------
		integrated <- RunUMAP(integrated, dims = 1:10)
		DimPlot(integrated, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.5)
		head(integrated@meta.data)
		integrated
		"

		#Annotation Diagnostics----------------------------------------------------------------------------
		pred.grun <- SingleR(test=as.SingleCellExperiment(integrated), ref=ref, assay.type.test=1, labels=ref$label.main, de.method="wilcox")
		table(pred.grun$labels)
		# plotScoreHeatmap(pred.grun)
		# plotDeltaDistribution(pred.grun, ncol = 3)

		# all.markers <- metadata(pred.grun)$de.genes
		# sceG$labels <- pred.grun$labels

		# Next, switch the identity class of all cells to reflect replicate ID
		Idents(integrated) <- "SingleR.calls"

		# How many cells are in each cluster
		table(Idents(integrated))

		# How many cells are in each cluster
		integrated_counts <- table(Idents(integrated))
		head(integrated_counts)

		integrated_long <- as.data.frame(integrated_counts)
		head(integrated_long)
		colnames(integrated_long) <- c("cell_type", "Cell_count")
		head(integrated_long)
		rownames(integrated_long) <- integrated_long$x
		integrated_long

		# Create ggplot2 plot scaled to 1.00
		bar_counts <- ggplot(integrated_long, aes(x = cell_type, y = Cell_count, fill = cell_type)) + 
			                    geom_bar(stat = "identity", color ="black") + scale_x_discrete(name ="Cell type") + 
					                        theme(legend.position="right") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

							# Draw ggplot2 plot scaled to 1.00
							bar_counts

							# save annotated file
							saveRDS(integrated, file = "annotation_integrated.rds")

