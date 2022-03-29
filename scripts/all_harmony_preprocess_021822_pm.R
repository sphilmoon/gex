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
library(harmony)
#############################################
# 01 PRE-PROCESSING
#############################################
# parameters for pre-processing are flexible.

# toc: Table of counts (filtered). Just those columns of \code{tod} that contain cells.
# tod: Table of droplets (raw).  A matrix with columns being each droplet and rows each gene.

# REMOVE BACKGROUND IN adj_4dpi
adj_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj/filtered_feature_bc_matrix")
adj_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj/raw_feature_bc_matrix")
adj_sc = SoupChannel(adj_tod, adj_toc)
adj_sc = SoupChannel(adj_tod, adj_toc, calcSoupProfile = FALSE)
adj_sc = estimateSoup(adj_sc)
adj_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj/clusters.csv")
adj_sc = setClusters(adj_sc, setNames(adj_metadata$Cluster, rownames(adj_metadata)))
# Estimate rho
adj_sc = autoEstCont(adj_sc)
# Clean the data
adj_out = adjustCounts(adj_sc)

# REMOVE BACKGROUND IN pyk_4dpi
pyk_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk/filtered_feature_bc_matrix")
pyk_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk/raw_feature_bc_matrix")
pyk_sc = SoupChannel(pyk_tod, pyk_toc)
pyk_sc = SoupChannel(pyk_tod, pyk_toc, calcSoupProfile = FALSE)
pyk_sc = estimateSoup(pyk_sc)
pyk_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk/clusters.csv")
pyk_sc = setClusters(pyk_sc, setNames(pyk_metadata$Cluster, rownames(pyk_metadata)))
# Estimate rho
pyk_sc = autoEstCont(pyk_sc)
# Clean the data
pyk_out = adjustCounts(pyk_sc)

# REMOVE BACKGROUND IN naive
naive_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/Naive_outs/naive_filtered_feature_bc_matrix")
naive_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/Naive_outs/naive_raw_feature_bc_matrix")
naive_sc = SoupChannel(naive_tod, naive_toc)
naive_sc = SoupChannel(naive_tod, naive_toc, calcSoupProfile = FALSE)
naive_sc = estimateSoup(naive_sc)
naive_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/Naive_outs/analysis/clustering/graphclust/clusters.csv")
naive_sc = setClusters(naive_sc, setNames(naive_metadata$Cluster, rownames(naive_metadata)))
# Estimate rho
naive_sc = autoEstCont(naive_sc)
# Clean the data
naive_out = adjustCounts(naive_sc)

# REMOVE BACKGROUND IN adj_only
adj_only_soupx = load10X('home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj_only/')
adj_only_soupx = autoEstCont(adj_only_soupx)
adj_only_out = adjustCounts(adj_only_soupx)

# REMOVE BACKGROUND IN pyk_only
pyk_only_soupx = load10X('home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk_only/')
pyk_only_soupx = autoEstCont(pyk_only_soupx)
pyk_only_out = pykustCounts(pyk_only_soupx)


# Set up adj object
adj <- CreateSeuratObject(counts = adj_out, project = "adj", min.cells = 5)
adj$group <- "adj"
adj <- subset(adj, subset = nFeature_RNA > 100)
adj <- NormalizeData(adj, verbose = FALSE)
adj <- FindVariableFeatures(adj, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from adj
counts.adj <- GetAssayData(adj, assay = "RNA")
counts.adj <- counts.adj[-(which(rownames(counts.adj) %in% c("Gm42418", "AY036118"))),]
adj <- subset(adj, features = rownames(counts.adj))

# Set up pyk object
pyk <- CreateSeuratObject(counts = pyk_out, project = "pyk", min.cells = 5)
pyk$group <- "pyk"
pyk <- subset(pyk, subset = nFeature_RNA > 100)
pyk <- NormalizeData(pyk, verbose = FALSE)
pyk <- FindVariableFeatures(pyk, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from pyk
counts.pyk <- GetAssayData(pyk, assay = "RNA")
counts.pyk <- counts.pyk[-(which(rownames(counts.pyk) %in% c("Gm42418", "AY036118"))),]
pyk <- subset(pyk, features = rownames(counts.pyk))

# Set up naive object
naive <- CreateSeuratObject(counts = naive_out, project = "naive", min.cells = 5)
naive$group <- "naive"
naive <- subset(naive, subset = nFeature_RNA > 100)
naive <- NormalizeData(naive, verbose = FALSE)
naive <- FindVariableFeatures(naive, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from naive
counts.naive <- GetAssayData(naive, assay = "RNA")
counts.naive <- counts.naive[-(which(rownames(counts.naive) %in% c("Gm42418", "AY036118"))),]
naive <- subset(naive, features = rownames(counts.naive))

# Set up adj_only object
adj_only <- CreateSeuratObject(counts = adj_only_out, project = "adj_only", min.cells = 5)
adj_only$group <- "adj_only"
adj_only <- subset(adj_only, subset = nFeature_RNA > 100)
adj_only <- NormalizeData(adj_only, verbose = FALSE)
adj_only <- FindVariableFeatures(adj_only, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from adj_only
counts.adj_only <- GetAssayData(adj_only, assay = "RNA")
counts.adj_only <- counts.adj_only[-(which(rownames(counts.adj_only) %in% c("Gm42418", "AY036118"))),]
adj_only <- subset(adj_only, features = rownames(counts.adj_only))

# Set up pyk_only object
pyk_only <- CreateSeuratObject(counts = pyk_only_out, project = "pyk_only", min.cells = 5)
pyk_only$group <- "pyk_only"
pyk_only <- subset(pyk_only, subset = nFeature_RNA > 100)
pyk_only <- NormalizeData(pyk_only, verbose = FALSE)
pyk_only <- FindVariableFeatures(pyk_only, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from pyk_only
counts.pyk_only <- GetAssayData(pyk_only, assay = "RNA")
counts.pyk_only <- counts.pyk_only[-(which(rownames(counts.pyk_only) %in% c("Gm42418", "AY036118"))),]
pyk_only <- subset(pyk_only, features = rownames(counts.pyk_only))





# QUALITY CONTROL adj
adj[["percent.MT"]] <- PercentageFeatureSet(adj, pattern = "^mt-")
VlnPlot(adj, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(adj, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(adj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 
adj <- subset(adj, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(adj, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL pyk
pyk[["percent.MT"]] <- PercentageFeatureSet(pyk, pattern = "^mt-")
VlnPlot(pyk, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(pyk, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(pyk, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
pyk <- subset(pyk, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(pyk, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL naive
naive[["percent.MT"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
naive <- subset(naive, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL adj_only
adj_only[["percent.MT"]] <- PercentageFeatureSet(adj_only, pattern = "^mt-")
VlnPlot(adj_only, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(adj_only, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(adj_only, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
adj_only <- subset(adj_only, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(adj_only, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL pyk_only
pyk_only[["percent.MT"]] <- PercentageFeatureSet(pyk_only, pattern = "^mt-")
VlnPlot(pyk_only, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(pyk_only, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(pyk_only, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
pyk_only <- subset(pyk_only, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(pyk_only, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)




### Remove batch effects by Harmony
### iterates until similar cells are clustered using PCA for dimensionality reduction.

combine <- merge(naive, y = list(adj, pyk), add.cell.ids = c("naive", "adj", "pyk", "adj_only", "pyk_only"), project = "PYK")
combine
combine <- NormalizeData(combine) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = combine, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = combine, features = "PC_2", group.by = "orig.ident",  pt.size = .1)
plot_grid(p1,p2)
options(repr.plot.height = 2.5, repr.plot.width = 6)
combine <- combine %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(combine, 'harmony')
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = combine, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = combine, features = "PC_1", group.by = "orig.ident",  pt.size = .1)
plot_grid(p1,p2)

ElbowPlot(combine, n=50)

combine <- combine %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(combine, reduction = "umap", group.by = "orig.ident", pt.size = .1)
##Check if cells cluster by cell cycles related genes
##If no, go to line 181-182 (save object)
##If yes, run from line 126 - 147  then run line 181-182 to save object
##Assign Cell-Cycle Scores

FeaturePlot(combine, features = c("Cd19", "Sdc1", "Ly6g", "Ly6c2", "Hmox1", "C1qa", "Itgax", "Klrb1c", "Cd4", "Cd8a", "Mki67", "Adgre1", "Aicda"))
DefaultAssay(combine) <- "RNA"
combine <- CellCycleScoring(combine, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(combine[[]])

# regressing out the difference between the G2M and S phase scores
combine$CC.Difference <- combine$S.Score - combine$G2M.Score
head(combine[[]])
combine <- ScaleData(combine, vars.to.regress = "CC.Difference", features = rownames(combine))

# after regressing, re-run clustering using harmony
combine <- combine %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)

# plot umap
DimPlot(combine, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(combine, reduction = "umap", group.by = "orig.ident")
DimPlot(combine, reduction = "umap", label = TRUE)

# SAVE
saveRDS(combine, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/all_harmony_preprocessing_021622_pm.rds")


