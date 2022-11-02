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
adj_4dpi_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj_4dpi/filtered_feature_bc_matrix")
adj_4dpi_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj_4dpi/raw_feature_bc_matrix")
adj_4dpi_sc = SoupChannel(adj_4dpi_tod, adj_4dpi_toc)
adj_4dpi_sc = SoupChannel(adj_4dpi_tod, adj_4dpi_toc, calcSoupProfile = FALSE)
adj_4dpi_sc = estimateSoup(adj_4dpi_sc)
adj_4dpi_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/adj_4dpi/clusters.csv")
adj_4dpi_sc = setClusters(adj_4dpi_sc, setNames(adj_4dpi_metadata$Cluster, rownames(adj_4dpi_metadata)))
# Estimate rho
adj_4dpi_sc = autoEstCont(adj_4dpi_sc)
# Clean the data
adj_4dpi_out = adjustCounts(adj_4dpi_sc)

# REMOVE BACKGROUND IN pyk_4dpi
pyk_4dpi_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk_4dpi/filtered_feature_bc_matrix")
pyk_4dpi_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk_4dpi/raw_feature_bc_matrix")
pyk_4dpi_sc = SoupChannel(pyk_4dpi_tod, pyk_4dpi_toc)
pyk_4dpi_sc = SoupChannel(pyk_4dpi_tod, pyk_4dpi_toc, calcSoupProfile = FALSE)
pyk_4dpi_sc = estimateSoup(pyk_4dpi_sc)
pyk_4dpi_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/dnalink_fastq/pyk_4dpi/clusters.csv")
pyk_4dpi_sc = setClusters(pyk_4dpi_sc, setNames(pyk_4dpi_metadata$Cluster, rownames(pyk_4dpi_metadata)))
# Estimate rho
pyk_4dpi_sc = autoEstCont(pyk_4dpi_sc)
# Clean the data
pyk_4dpi_out = adjustCounts(pyk_4dpi_sc)

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

# REMOVE BACKGROUND IN adj_naive
adj_naive_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_adj/outs/per_sample_outs/Navie_adj/count/count/sample_feature_bc_matrix")
adj_naive_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_adj/outs/multi/count/raw_feature_bc_matrix")
adj_naive_sc = SoupChannel(adj_naive_tod, adj_naive_toc)
adj_naive_sc = SoupChannel(adj_naive_tod, adj_naive_toc, calcSoupProfile = FALSE)
adj_naive_sc = estimateSoup(adj_naive_sc)
adj_naive_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_adj/outs/per_sample_outs/Navie_adj/count/count/analysis/clustering/graphclust/clusters.csv")
adj_naive_sc = setClusters(adj_naive_sc, setNames(adj_naive_metadata$Cluster, rownames(adj_naive_metadata)))
# Estimate rho
adj_naive_sc = autoEstCont(adj_naive_sc)
# Clean the data
adj_naive_out = adjustCounts(adj_naive_sc)

# REMOVE BACKGROUND IN pyk_naive
pyk_naive_toc = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_PYK/
                        outs/per_sample_outs/Navie_PYK/count/sample_feature_bc_matrix")
pyk_naive_tod = Read10X(data.dir ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_PYK/
                        outs/multi/count/raw_feature_bc_matrix")
pyk_naive_sc = SoupChannel(pyk_naive_tod, pyk_naive_toc)
pyk_naive_sc = SoupChannel(pyk_naive_tod, pyk_naive_toc, calcSoupProfile = FALSE)
pyk_naive_sc = estimateSoup(pyk_naive_sc)
pyk_naive_metadata <-read.csv(file ="/home/bmrc/Public/phil_ubuntu/sc/github/gex_jupyter/raw_data/pyk_datasets/adj_pyk_naive/Navie_PYK
                              /outs/per_sample_outs/Navie_PYK/count/analysis/clustering/graphclust/clusters.csv")
pyk_naive_sc = setClusters(pyk_naive_sc, setNames(pyk_naive_metadata$Cluster, rownames(pyk_naive_metadata)))
# Estimate rho
pyk_naive_sc = autoEstCont(pyk_naive_sc)
# Clean the data
pyk_naive_out = adjustCounts(pyk_naive_sc)

### using both count matrix and clustering results for my 10X dataset
# Set up adj_4dpi object
adj_4dpi <- CreateSeuratObject(counts = adj_4dpi_out, project = "adj_4dpi", min.cells = 5)
adj_4dpi$group <- "adj_4dpi"
adj_4dpi <- subset(adj_4dpi, subset = nFeature_RNA > 100)
adj_4dpi <- NormalizeData(adj_4dpi, verbose = FALSE)
adj_4dpi <- FindVariableFeatures(adj_4dpi, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from adj_4dpi
counts.adj_4dpi <- GetAssayData(adj_4dpi, assay = "RNA")
counts.adj_4dpi <- counts.adj_4dpi[-(which(rownames(counts.adj_4dpi) %in% c("Gm42418", "AY036118"))),]
adj_4dpi <- subset(adj_4dpi, features = rownames(counts.adj_4dpi))

# Set up pyk_4dpi object
pyk_4dpi <- CreateSeuratObject(counts = pyk_4dpi_out, project = "pyk_4dpi", min.cells = 5)
pyk_4dpi$group <- "pyk_4dpi"
pyk_4dpi <- subset(pyk_4dpi, subset = nFeature_RNA > 100)
pyk_4dpi <- NormalizeData(pyk_4dpi, verbose = FALSE)
pyk_4dpi <- FindVariableFeatures(pyk_4dpi, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from pyk_4dpi
counts.pyk_4dpi <- GetAssayData(pyk_4dpi, assay = "RNA")
counts.pyk_4dpi <- counts.pyk_4dpi[-(which(rownames(counts.pyk_4dpi) %in% c("Gm42418", "AY036118"))),]
pyk_4dpi <- subset(pyk_4dpi, features = rownames(counts.pyk_4dpi))

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
adj_naive <- CreateSeuratObject(counts = adj_naive_out, project = "adj_naive", min.cells = 5)
adj_naive$group <- "adj_naive"
adj_naive <- subset(adj_naive, subset = nFeature_RNA > 100)
adj_naive <- NormalizeData(adj_naive, verbose = FALSE)
adj_naive <- FindVariableFeatures(adj_naive, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from adj_naive
counts.adj_naive <- GetAssayData(adj_naive, assay = "RNA")
counts.adj_naive <- counts.adj_naive[-(which(rownames(counts.adj_naive) %in% c("Gm42418", "AY036118"))),]
adj_naive <- subset(adj_naive, features = rownames(counts.adj_naive))

# Set up pyk_naive object
pyk_naive <- CreateSeuratObject(counts = pyk_naive_out, project = "pyk_naive", min.cells = 5)
pyk_naive$group <- "pyk_naive"
pyk_naive <- subset(pyk_naive, subset = nFeature_RNA > 100)
pyk_naive <- NormalizeData(pyk_naive, verbose = FALSE)
pyk_naive <- FindVariableFeatures(pyk_naive, selection.method = "vst", nfeatures = 2000)
## Remove Gm42418 (rRNA gene) and AY036118 (long non-coding RNA) from pyk_naive
counts.pyk_naive <- GetAssayData(pyk_naive, assay = "RNA")
counts.pyk_naive <- counts.pyk_naive[-(which(rownames(counts.pyk_naive) %in% c("Gm42418", "AY036118"))),]
pyk_naive <- subset(pyk_naive, features = rownames(counts.pyk_naive))



# QUALITY CONTROL adj_4dpi
adj_4dpi[["percent.MT"]] <- PercentageFeatureSet(adj_4dpi, pattern = "^mt-")
VlnPlot(adj_4dpi, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot1 <- FeatureScatter(adj_4dpi, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot2 <- FeatureScatter(adj_4dpi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 
adj_4dpi <- subset(adj_4dpi, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(adj_4dpi, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL pyk_4dpi
pyk_4dpi[["percent.MT"]] <- PercentageFeatureSet(pyk_4dpi, pattern = "^mt-")
VlnPlot(pyk_4dpi, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(pyk_4dpi, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(pyk_4dpi, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
pyk_4dpi <- subset(pyk_4dpi, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(pyk_4dpi, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL naive
naive[["percent.MT"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
naive <- subset(naive, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL adj_naive
adj_naive[["percent.MT"]] <- PercentageFeatureSet(adj_naive, pattern = "^mt-")
VlnPlot(adj_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(adj_naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(adj_naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
adj_naive <- subset(adj_naive, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(adj_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

# QUALITY CONTROL pyk_naive
pyk_naive[["percent.MT"]] <- PercentageFeatureSet(pyk_naive, pattern = "^mt-")
VlnPlot(pyk_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)
plot3 <- FeatureScatter(pyk_naive, feature1 = "nCount_RNA", feature2 = "percent.MT")
plot4 <- FeatureScatter(pyk_naive, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4 
pyk_naive <- subset(pyk_naive, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.MT < 10)
VlnPlot(pyk_naive, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)


### Remove batch effects by Harmony between two different batch time for sequencing
### iterates until similar cells are clustered using PCA for dimensionality reduction.

all5 <- merge(naive, y = list(adj_4dpi, pyk_4dpi, adj_naive, pyk_naive), add.cell.ids = c("naive", "adj_4dpi", "pyk_4dpi", "adj_naive", "pyk_naive"), 
              project = "PYK_vaccine")
all5
all5 <- NormalizeData(all5) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = all5, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = all5, features = "PC_2", group.by = "orig.ident",  pt.size = .1)
plot_grid(p1,p2)
options(repr.plot.height = 2.5, repr.plot.width = 6)
all5 <- all5 %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(all5, 'harmony')
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = all5, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = all5, features = "PC_1", group.by = "orig.ident",  pt.size = .1)
plot_grid(p1,p2)

ElbowPlot(all5, n=50)

all5 <- all5 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 1) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)
DimPlot(all5, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(all5, reduction = "umap", split.by = "orig.ident", label = FALSE)

##Check if cells cluster by cell cycles related genes
##If no, go to line 181-182 (save object)
##If yes, run from line 126 - 147  then run line 181-182 to save object
##Assign Cell-Cycle Scores

FeaturePlot(all5, features = c("Cd19", "Sdc1", "Ly6g", "Ly6c2", "Hmox1", "C1qa", "Itgax", "Klrb1c", 
                               "Cd4", "Cd8a", "Mki67", "Adgre1", "Aicda","Ear2", "Ace", "S100a8",
                               "Cd3e", "Foxp3", "Ptprc"))
DefaultAssay(all5) <- "RNA"
all5 <- CellCycleScoring(all5, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(all5[[]])

# regressing out the difference between the G2M and S phase scores
all5$CC.Difference <- all5$S.Score - all5$G2M.Score
head(all5[[]])
all5 <- ScaleData(all5, vars.to.regress = "CC.Difference", features = rownames(all5))

# resolution at 0.5 -> 6 Bsubset clutsers
# resolution at 0.75 -> 9 Bsubset clutsers
# resolution at 0.1 -> 11 Bsubset clutsers


# Everytime when perform normalization, always re-cluster. 
# after regressing, re-run clustering using UMAP with harmony
ElbowPlot(all5, n=50)
all5 <- all5 %>% 
  RunUMAP(reduction = "harmony", dims = 1:40) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:40) %>% 
  FindClusters(resolution = 0.75) %>% 
  identity()
options(repr.plot.height = 4, repr.plot.width = 10)

# plot umap
DimPlot(all5, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(all5, reduction = "umap", group.by = "orig.ident")
DimPlot(all5, reduction = "umap", label = TRUE)


# Load ImmGene reference
ref <- ImmGenData()
ref

# Assigning general annotation
annotation <- SingleR(test=as.SingleCellExperiment(all5), ref=ref, assay.type.test=1, 
                      labels=ref$label.main)
annotation
head(all5@meta.data)
table(annotation$labels)

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
all5$SingleR.pruned.calls <- annotation$pruned.labels
all5$SingleR.calls <- annotation$labels

# Run UMAP (as well as pca and harmony)
all5 <- RunUMAP(all5, dims = 1:28)
DimPlot(all5, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.1)
head(all5@meta.data)
all5

# Annotation Diagnostics
pred.grun <- SingleR(test=as.SingleCellExperiment(all5), ref=ref, assay.type.test=1, 
                     labels=ref$label.main, de.method="wilcox")
table(pred.grun$labels)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(all5) <- "SingleR.calls"

# How many cells are in each cluster
all5_counts <- table(Idents(all5))
head(all5_counts, n =20)

# remove unnecessary clusters
all5_removal <- subset((all5), 
                       idents = c("ILC", "Tgd", "Stem cells", "B cells, pro",
                                  "Basophils", "Eosinophils",
                                  "Microglia", "Epithelial cells", "Mast cells",
                                  "Endothelial cells"), 
                       invert = T)
# now visualize the heterogenity
DimPlot(all5_removal, reduction = "umap", label = T, label.size = 4.5) + NoLegend()

all5_removal_counts <- table(Idents(all5_removal), all5_removal$group)
head(all5_removal_counts, n =10)
write.csv(all5_removal_counts, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/all5_removal_counts_081522_PM.csv")


# general markers for the entire heterogenity 
# 1. macrophages Nrp1, C1qa, Ms4a7
# 2. t cells Cd4, Cd8a, Cd3e, Foxp3
# 3. B cells Cd19, Cd79a, Cd79b, Pax5
# 4. NK Ncr1, Fcgr3a
# 5. DC Ccr7, Tifab
# 6. NKT Klrb1c, Trdc, Cd3d, Ccl5
# 7. Neutrophils Ly6c2, Ear2
# 8. Monocytes Cd14, Cd16, S100a9, Lyz


hetero_genes <- c("Cd19", "Cd79a", "Cd79b","Pax5",
                  "Cd4", "Cd8a",  "Cd3d", "Il7r", "Cd3e", "Trac", "Xcl1", 
                  "Ncr1", "Klrb1c", "Klrk1", "Klrg1", "Gzma", "Gzmb", "Sell", "Klrd1",
                  "Ly6g", "Cd14", "Csf3r", "Fcgr3", "Cd177",
                  "Ccr7", "Clec9a", "Cst3", "Itgax",
                  "Ly6c2", "Nrp1", "Fcgr1",
                  "C1qa", "C1qb", "C1qc", "Ms4a7", "Adgre1")

# reorder clusters
all5_removal@active.ident <- factor(all5_removal@active.ident, 
                        levels=c("Macrophages", 
                                 "Monocytes",
                                 "DC", 
                                 "Neutrophils", 
                                 "NK cells", 
                                 "NKT", 
                                 "T cells", 
                                 "B cells"))

DotPlot(object = all5_removal, features = hetero_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# SAVE 
# saved
saveRDS(all5, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/01_all5_res75_dr40_prepro_harmony_080822_pm.rds")
# read annotated file
all5 <- readRDS(file= "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/01_all5_res75_dr40_prepro_harmony_080822_pm.rds")

#### 02_B cell subsetting using Immgen

# Next, switch the identity class of all cells to reflect replicate ID
# now using ImmGen identification label
Idents(all5) <- "SingleR.calls"

# How many cells are in each cluster
all5_counts <- table(Idents(all5))
head(all5_counts, n =20)

# plot umap
DimPlot(all5, reduction = "umap", split.by = "orig.ident", label = FALSE)
DimPlot(all5, reduction = "umap", group.by = "orig.ident")
DimPlot(all5, reduction = "umap", label = TRUE)
# using SingleR annotation
DimPlot(all5, reduction = "umap", group.by = "SingleR.calls", split.by = "orig.ident", label = TRUE, pt.size = 0.1)
DimPlot(all5, reduction = "umap", group.by = "SingleR.calls", label = T, pt.size = 0.1)

# only 63 pro B cell were annotated from Immgen database, so it is excluded.
bcellsubset_immgen <- subset((all5), idents = c("B cells"))

# B cells
bcell <- c("Cd19", "Ms4a1", "Cd79a", "Cd79b", "Ighd", "Ighm", "Jchain", "Cd1d1")
FeaturePlot(all5, features = bcell, label = T)
DotPlot(all4, features = bcell, cols = c("blue", "red"))

?subset

# How many cells are in each cluster
table(Idents(bcellsubset_immgen))
bcellsubset_immgen
head(bcellsubset_immgen@meta.data)

# Visualize and just select B cells before subsetting.
DimPlot(bcellsubset_immgen, reduction = "umap",  label = T, pt.size = 1)
DimPlot(bcellsubset_immgen, reduction = "umap", group.by = "orig.ident")

## Neccessary for B cell subset clutering
################**************################ 
# Now normalize B cells and cluster UMAP for subsetting.
bcellsubset_immgen <- NormalizeData(bcellsubset_immgen) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA(verbose = FALSE)
ElbowPlot(bcellsubset_immgen, n=50)
bcellsubset_immgen <- bcellsubset_immgen %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
bcellsubset_immgen <- bcellsubset_immgen %>% 
  RunUMAP(reduction = "harmony", dims = 1:32) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:32) %>% 
  FindClusters(resolution = 0.35) %>% 
  identity()
# Visualize the newly normalized and clustered B cells
DimPlot(bcellsubset_immgen, reduction = "umap",  label = T, pt.size = 1)
################**************################

# save annotated file 
# saved
saveRDS(bcellsubset_immgen, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/final_all5_bcellsubset_080822_pm.rds")
bcellsubset_immgen <- readRDS(file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/final_all5_bcellsubset_080822_pm.rds")

# general annotation using B cell canonical markers
FeaturePlot(bcellsubset_immgen, features = bcell, label = T)
DotPlot(bcellsubset_immgen, features = bcell, cols = c("blue", "red"))

# filter out non-B cells

dc <- c("Itgam", "Itgax", "Il12a", "Ifng", "Trbc1", "Ccr7", "Fcer1a", "Cst3")
FeaturePlot(bcellsubset_immgen, features = dc, label = T)
DotPlot(bcellsubset_immgen, features = dc, cols = c("blue", "red"))

tcell <- c("Trac", "Trbc1", "Cd3e", "Cd4", "Cd8a", "Foxp3") 
FeaturePlot(bcellsubset_immgen, features = tcell, label = T) # cluster 9

# Macrophages
FeaturePlot(bcellsubset_immgen, features = c("C1qa", "Ms4a7")) 
# NK cells
FeaturePlot(bcellsubset_immgen, features = c("Ncr1", "Klrb1c", "Nkg7"))
# Neutrophils/Monocytes
FeaturePlot(bcellsubset_immgen, features = c("Ly6c2", "Ear2", "Fcgr3a", "Ms4a7"))


# visualize B subsets
DimPlot(bcellsubset_immgen, reduction = "umap", label = T)


# MBC ?
mbc_genes <- c("Cr2", "Cd52", "Vim", "Parm1", "Cd55", "Foxp1", "Cxcr4", # Riedel_nature_2020
               "Cxcr3", "Ptpn22", "Cd9", "S1pr1", "Il10", "Il10ra",
               "Ighg1", "Ighg2b", "Ighg2c", "Ighg3",
               "Pdcd1", "Pdcd1lg2")
FeaturePlot(bcellsubset_immgen, features = mbc_genes, label = T)

DotPlot(object = bcellsubset_immgen, features = mbc_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")



# remove non-B cells or multiplets
bcellsubset_removal <- subset((bcellsubset_immgen), 
                              idents = c(6, 7, 8, 9, 10, 12, 13), 
                              invert = T)
DimPlot(bcellsubset_removal, label = T, pt.size = 0.1)


### Annotation ###

# cluster 0
# Follicular B cells
# Ighd, Fcer2a, (CD23, IgE), Cd55

fob_genes <- c("Ighd", "Fcer2a", "Cd19", "Cr2", "Cd1d1", "Ms4a1", "Cd40",
               "Cxcr5", "Sell", "Cd93", "H2-Aa", "Shisa5", "H2-Eb1", "H2-Ab1",
               "Fcmr", "Cd55", "Fchsd2", "Irf1")
DotPlot(object = bcellsubset_removal, features = fob_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

FeaturePlot(bcellsubset_removal, features = fob_genes, label = T)
VlnPlot(bcellsubset_removal, features = fob_genes, pt.size = 0.3, 
        group.by = "orig.ident")

# cluster 1
# T1 B cells
# Cd93, Cd24a, Ly6d, Ebf1, Ms4a1, Vpreb3 high expression
# Cr2 (Cd21), Fcer2a (Cd23) absence

t1b <- c("Cd19", "Ighm", "Cd93", "Ighd", "Cr2", "Fcer2a", "Cd24a",
         "Ly6d", "Ebf1", "Ms4a1", "Vpreb3", "Iglc1", "Klf2")

DotPlot(object = bcellsubset_removal, features = t1b) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")
FeaturePlot(bcellsubset_removal, features = t1b, label = T)

c2_genes <- c("Ly6d", "Ms4a1", "Ighd", "Ighm", "Fcer2a", "Cd19", "Cr2", "Cd93",
              "Cd24a")

# cluster 2
# T2/T3
t2_t3 <- c("Ly6d", "Ms4a1", "Ighd", "Ighm", "Fcer2a", "Cd19", "Cr2", "Cd93",
           "Cd24a")
DotPlot(object = bcellsubset_removal, features = t2_t3) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# cluster 3
# Marginal Zone B cells
# Cr2 and Cd1d1 major genes

mzb <- c("Ly6d", "Ebf1", "Cr2", "Ms4a1", "Cd19", "Cd1d1", "Ccr2", "Zbtb32", 
              "Rsu1", "Dph5", "Zfp36l1")
DotPlot(object = b_cells_subset, features = mzb) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# cluster 4
# B reg
# Apoe, Cd74, Bach2
breg <- c("Fcrl5", "Zbtb20", "Ccdc28b", "Cd9", "Ptpn22", "Atf3", "Il10", "Apoe")
FeaturePlot(bcellsubset_removal, features = breg, label = T)
DotPlot(object = bcellsubset_removal, features = breg) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# B1 B cells
# Zbtb32, Ahnak, Vim, S100a6, Bhlhe41
b1b_genes <- c("Itgam", "Itgax", "Tbx21", "S100a6", "Bhlhe41", "Fcrl5",
               "Cdk1", "Ahnak", "Vim")
FeaturePlot(bcellsubset_removal, features = b1b_genes, label = T)
DotPlot(object = bcellsubset_removal, features = b1b_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")


# prepro
prepro <- c("Ebf1", "Ybx3", "Cd93", "Il2ra", "Spn", "Cd43", "Cd25", "Il7r", "Flt3")
FeaturePlot(bcellsubset_removal, features = prepro, label = T)
DotPlot(object = bcellsubset_removal, features = prepro) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# cluster 5
# GC
# Aicda, Mki67, Cdk1, Cdc20, Pcna, Pclaf
gc <- c("Aicda", "Mki67", "Cdk1", "Cdc20", "Pcna", "Pclaf", "Bcl6", "Acta2", "Col1a1")
FeaturePlot(bcell_subset_final, features = gc, label = T)
DotPlot(bcellsubset_removal, features = gc, cols = c("blue", "red"))
DimPlot(bcell_subset_final, features = gc, reduction = "umap", group.by = "orig.ident") + DarkTheme()

cluster_gc <- FindMarkers(bcell_subset_final, ident.1 = "GC", only.pos = T, min.pct = 0.25)
head(cluster_gc, n = 40)


# cluster 11
# Plasma B cells
# Sdc1, Xbp1, Jchain, Ighm high expression
pc <- c("Sdc1", "Xbp1", "Ighm", "Jchain", "Ighg1", "Ighg2b", "Ighg3",
        "Igha", "Prdm1")
FeaturePlot(bcellsubset_removal, features = pc, label = T)
DotPlot(bcellsubset_removal, features = pc, cols = c("blue", "red"))

# saved
saveRDS(bcellsubset_removal, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/03_bcellsubset_removal_081522_pm.rds")

# x <- readRDS(file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/03_bcellsubset_removal_081522_pm.rds")

########################
### FINAL ANNOTATION ###
########################
# Naming clusters based on canonical markers
new_bcluster <- c("FoB", "T1 B", "T2/T3 B", "MzB", "B1 B", "GC", "6", 
                  "7_ifn", "8_dc", "9_macro", "10", "PC", "12_dc", "13_dc")
names(new_bcluster) <- levels(bcellsubset_immgen)
bcellsubset_immgen <- RenameIdents(bcellsubset_immgen, new_bcluster)

# remove non-B cells or multiplets
bcell_subset_final <- subset((bcellsubset_immgen), 
                            idents = c("6", "7_ifn", "8_dc", "9_macro", "10", 
                                       "12_dc", "13_dc"), 
                            invert = T)

DimPlot(bcell_subset_final, label = T, pt.size = 0.1)
DimPlot(bcell_subset_final, reduction = "umap", group.by = "orig.ident") + DarkTheme()
DimPlot(bcell_subset_final, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.1) 

#### 04 ####
# saved
saveRDS(bcell_subset_final, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/04_bcell_subset_final_081522_pm.rds")
readRDS(bcell_subset_final, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/04_bcell_subset_final_081522_pm.rds")
# Daliya's suggestion
DimPlot(bcell_subset_final[,WhichCells(bcell_subset_final, expression=orig.ident != "naive")], reduction = "umap", split.by = "orig.ident", label = T, pt.size = 0.5)


# B cell subset dot plot
bcell_final_dotplot <- c("Cd19", "Cd79a", "Cd79b", "Pax5", "Ly6d",
                         "Ighd", "Cd55", "Fcer2a",
                         "Cd93", "Cd24a", "Ms4a1", "Vpreb3", 
                         "Cr2", "Cd1d1", "Bcl6",
                         "Ebf1", "Apoe", "Plac8", "Zbtb20", "Zbtb32", "Ahnak", "Vim", "S100a6", "Fcrl5", "Ly86", "Bhlhe41",
                         "Aicda", "Mki67", "Cdk1", "Cdc20", "Pcna", "Pclaf", "Ube2c", 
                         "Sdc1", "Xbp1", "Prdm1", "Irf4", "Jchain", "Iglc1", "Ighm", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3")
DotPlot(object = bcell_subset_final, features = bcell_final_dotplot) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")


# new mbc dot plot
new_mbc_genes <- c("Id2", "S100a4", "Pdcd1", "Pdcd1lg2", "Parm1", # Riedel_nature_2020
                   "Cd28", "Cd40", "Cd44", "Cd247", "Zap70", "Arl4c",
                   "Itgam", "Itgax", "Tbx21", #atMBC
                   "Cd19", "Cd38", "Sdc1") 

gc_genes <- c("Aicda", "Mki67", "Cdk1", "Cdc20", "Pcna", "Pclaf", "Ube2c", "Bcl6")
pc_genes <- c("Sdc1", "Xbp1", "Prdm1", "Irf4", "Jchain", "Iglc1", "Ighm", 
              "Ighg1", "Ighg2b", "Ighg2c", "Ighg3")


# figures
DotPlot(object = bcell_subset_final, features = new_mbc_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# mbc feature plot
FeaturePlot(bcell_subset_final, features = new_mbc_genes, label = F) + scale_fill_gradientn(colors = c("blue", "white", "red"))
# gc feature plot
FeaturePlot(bcell_subset_final, features = gc_genes, label = F) + scale_fill_gradientn(colors = c("blue", "white", "red"))
# pc feature plot
FeaturePlot(bcell_subset_final, features = pc_genes, label = F) + scale_fill_gradientn(colors = c("blue", "white", "red"))


dge_bcell_subset_final <- FindAllMarkers(bcell_subset_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.05)
top10_dge_bcell_subset_final <- dge_bcell_subset_final %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

DoHeatmap(bcell_subset_final, features = top10_dge_bcell_subset_final$gene, group.bar = TRUE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
DoHeatmap(bcell_subset_final, features = top10_dge_bcell_subset_final$gene, group.bar = TRUE) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

# subsetting GC only
gc_final <- subset((bcell_subset_final), 
                             idents = c("FoB", "T1 B", "T2/T3 B", "MzB", "B1 B", 
                                        "PC"), 
                             invert = T)
# subsetting PC only
pc_final <- subset((bcell_subset_final), 
                   idents = c("FoB", "T1 B", "T2/T3 B", "MzB", "B1 B", 
                              "GC"), 
                   invert = T)


# GO vln plots

# SHM
shm_genes <- c("Aicda", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3")

VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = shm_genes, pt.size = 0.2, split.by = "orig.ident",
        stack = T, flip=T, same.y.lims = T) 
        
VlnPlot(gc_final[,WhichCells(gc_final, expression=orig.ident != "naive")], 
        features = shm_genes, pt.size = 0.3, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

DoHeatmap(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = shm_genes) + scale_fill_gradientn(colors = c("blue", "white", "red"))

# Cell diff

cell_diff_gc_pyk_adj <- c("Mef2b", "Ppp4r2", "Mybl1", "Bcl6", "Gadd45b")

VlnPlot(gc_final[,WhichCells(gc_final, expression=orig.ident != "naive")], 
        features = cell_diff_gc_pyk_adj, pt.size = 0.2, split.by = "orig.ident",
        stack = T, flip=T, same.y.lims = T) 

# gene counts
GetAssayData(gc_final, slot="counts")[c("Aicda", "Bcl6"),]

cell_death_genes <- c("Bad", "Bax", "Bak1", "Puma", "Diablo", "Noxa", "Bim",
                      "Apaf", "Cyb561d2", "Endog", "Gadd45a", "Casp1",
                      "Casp2", "Casp3", "Casp8", 
                      "Casp11", "Casp14", "Traad", "Fas",
                      "Tnfsf10", "Bca3", "Bcl2l11", 
                      "Ripk1", "Ripk3", "Zbpi", "Nrg3")

cell_proliferation_genes_1 <- c("Mcl1", "Bcl2", "Traf3", 
                              "Tnfrsf17", "Tnfrsf13b", "Bcl2a1b", "Bcl10",
                              "Xiap", "Birc5", "Cflar", "Traf1", "Mki67",
                              "Anp32b", "Pmf1", "Stmn1", "Ccne2", 
                              "Ranbp1", "Cks1b", "Cenpf", "Cdc20", "Smc4")

cell_proliferation_genes_2 <- c("H2afx", "Cks2", "Top2a", "Smc2", "Spc24",
                              "Tuba1c", "Cdkn2d", "Cdca8", "Cdk2ap2", 
                              "Cdk1", "Cdk4", "Cdk9", "Cdk11b", "Cdk5rap3",
                              "Cdk13", "Cdk12")

cell_proliferation_genes_3 <- c("Ccni", "Ccnd3", "Ccnl2",
                              "Ccng2", "Ccnd2", "Ccnk", "Ccny", "Ccna2", 
                              "Ccnt1", "Ccnb2", "Ccnh", "Ccnb1", "Ccne1", 
                              "Cd40", "Tnfrsf13c")

mzb_fob_gc_pc_final <- subset((bcell_subset_final), 
                   idents = c("T1 B", "T2/T3 B", "B1 B"), invert = T)

# Aicda
VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = gc_genes, pt.size = 0.1, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

# Cell death
VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = cell_death_genes, pt.size = 0.1, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

DoHeatmap(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
          features = cell_death_genes) + scale_fill_gradientn(colors = c("blue", "white", "red"))

# Cell proliferation_1
VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = cell_proliferation_genes_1, pt.size = 0.1, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

# Cell proliferation_2
VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = cell_proliferation_genes_2, pt.size = 0.1, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

# Cell proliferation_3
VlnPlot(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
        features = cell_proliferation_genes_3, pt.size = 0.1, split.by = "orig.ident",
        stack=T, flip=T, same.y.lims = T) & 
  geom_boxplot(width = 0.2,position = position_dodge(0.9))

DoHeatmap(mzb_fob_gc_pc_final[,WhichCells(mzb_fob_gc_pc_final, expression=orig.ident != "naive")], 
          features = cell_proliferation_genes) + scale_fill_gradientn(colors = c("blue", "white", "red"))

############################
### Cell number analysis ###
############################
# How many cells are in each cluster
bcell_subset_final_count <- table(Idents(bcell_subset_final), bcell_subset_final$group)
head(bcell_subset_final_count, n =15)
# write.csv(bcell_subset_final_count, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/all5_counts_bcell_subset_final_count_081522_PM.csv")

# making bar graphs
bcell_subset_final_count <- table(Idents(bcell_subset_final), bcell_subset_final$group)
pt <- as.data.frame(bcell_subset_final_count)
pt$Var1 <- as.character(pt$Var1)
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "black")+
  xlab("Condition") +
  ylab("Cell number") +
  geom_text(aes(label = Freq), fontface = "bold", size =10, 
            position = position_stack(vjust = 0.5)) + theme_bw(base_size = 0) +
            theme(legend.position = c(0.5, 0.8), 
                  legend.title = element_text(face="bold",size=30),
                  legend.text = element_text(face="bold", size=25), 
                  legend.background = element_rect(size=0.8, color="black"),
                  axis.text.x= element_text(face="bold",size=25))

bcell_subset_final$group <- factor(bcell_subset_final$group, 
                                   levels = c("adj_naive", "pyk_naive", 
                                              "adj_4dpi", "pyk_4dpi", "naive"))

ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "black")+
  xlab("Condition") +
  ylab("Cell number") +
  geom_text(aes(label = Freq), fontface = "bold", size =10, 
            position = position_stack(vjust = 0.5)) + theme_bw(base_size = 0) +
  theme(legend.position = c(0.5, 0.8), 
        legend.title = element_text(face="bold",size=30),
        legend.text = element_text(face="bold", size=25), 
        legend.background = element_rect(size=0.8, color="black"),
        axis.text.x= element_text(face="bold",size=25))


# making pie charts
colnames(pt) <- c("Celltype", "Condition", "Frequency")
pt2 <- pt %>%  
  group_by(Condition) %>%
  add_tally(Frequency, name = "total") %>%
  ungroup() %>%
  mutate(label = sprintf("total = %d", total))

pt2 <- pt2 %>% mutate(total2 = paste0(case_when(Condition=='pyk_naive' ~ round(Frequency*83030, digits= 2),
                                                Condition=='pyk_4dpi' ~ round(Frequency*103770, digits= 2),
                                                Condition=="adj_4dpi" ~ round(Frequency*86770, digits=2),
                                                Condition=='adj_naive' ~ round(Frequency*81670, digits= 2))))
pt2 <- pt2 %>%  
  group_by(Condition) %>%
  add_tally(as.numeric(total2), name = "total3") %>%
  ungroup() %>%
  mutate(label2 = sprintf("total3 = %d", total3))

pt2 %>% ggplot(aes(x = total3/2, y = Frequency, fill = `Celltype`, width = total3)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  facet_wrap(~ Condition, strip.position = "bottom") +
  coord_polar("y", start = 0, direction = -1) +
  theme_bw(base_size = 0) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text = element_text(face = "bold", size = 25)) + NoLegend()


# 02
colnames(pt) <- c("Celltype", "Condition", "Frequency")
pt$Condition <- factor(pt$`Condition`, level = c("adj_naive", "adj_4dpi"))
ggplot(pt, aes(x = Condition, y = Frequency, fill = "Celltype")) +
  geom_bar(stat = "identity", color = 'black')+
  xlab("Condition") + ylab("Cell number") +
  geom_text(aes(label = Frequency), fontface = "bold",size = 10, position = position_stack(vjust = 0.5)) +
  theme(legend.position = c(0.25,0.8), legend.title = element_text(face="bold",size=30),
        legend.text = element_text(face="bold",size=25), legend.background = element_rect(size=0.8,color="black"),
        axis.text.x= element_text(face="bold",size=25), plot.title = element_text(size=25)) +
  labs(title='(D) Barplot of DC(Ev14dpi)')



# save all4 final B cell subset annotation
saveRDS(bcell_subset_final, 
        file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/all5_final_annotation_081422_pm.rds")



#######################################################################
# 02_new DGE_081522_pm
# As a default, Seurat uses the non-parametric Wilcoxon rank sum test. 
#######################################################################


# read previously pre-processed, dimensional analysis, and clustered file
# bcell_subset_final_dge <- readRDS(file="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/all5_final_annotation_081422_pm.rds")


# 01. PYK vs ADJ
# A. PC
# pyk_naive vs adj_naive within the PC cluster
# first, looking at the pyk-vaccinated induced genes
pc_dge_genes_pyk_adj <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                    ident.1 = "pyk_naive", 
                                    ident.2 = "adj_naive", 
                                    subset.ident = "PC", min.pct = 0.25,
                                    only.pos = TRUE)
head(pc_dge_genes_pyk_adj, n = 20)
dim(pc_dge_genes_pyk_adj)

write.csv(pc_dge_genes_pyk_adj, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/pc_dge_genes_pyk_adj_081522_pm.csv')
EnhancedVolcano(pc_dge_genes_pyk_adj,
                lab = rownames(pc_dge_genes_pyk_adj),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5,
                xlim = c(-3, 3),
                ylim = c(-0.1, 5),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'none',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')


# B. GC
# pyk_naive vs adj_naive within the GC cluster
# first, looking at the pyk-vaccinated induced genes
gc_dge_genes_pyk_adj <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                    ident.1 = "pyk_naive", 
                                    ident.2 = "adj_naive", 
                                    subset.ident = "GC", min.pct = 0.25,
                                    only.pos = F)
head(gc_dge_genes_pyk_adj, n = 20)
dim(gc_dge_genes_pyk_adj)

write.csv(gc_dge_genes_pyk_adj, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk_adj_081522_pm.csv')
EnhancedVolcano(gc_dge_genes_pyk_adj,
                lab = rownames(gc_dge_genes_pyk_adj),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c("Polr2g", "Ppp4r2", "H1f10", "Bcl6", "Nuggc", 
                "Aicda", "Pold4", "Cdk2ap2", "Ywhah", "Dynll1", "Chmp2a",
                "Bcl6", "Ppp4r2", "H1f10", "Mbd4", "Top1", "Elob", "Atp5g3",
                "Ifi203", "Foxp1", "Creg1", "Rpl4", "Eif3h", "Rbm3", "Ptprcap"),
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5,
                xlim = c(-3, 3),
                ylim = c(-0.1, 7),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

#selectLab = c('Atp5g3', 'Mef2b', 'Dynll1', 'Chmp2a',
 #             'Oaz1', 'Polr2g', 'Ppp4r2', 'Elob', 'Nuggc',
  #            'Mybl1', 'Top1', 'Gadd45b', 'Traf1', 'Uqcrc2',
   #           'H1fx', 'Ywhah', 'Aicda', 'Bcl6',
    #          "mt-Nd3", "Ifi203", "Foxp1", "Lars2"),


# 02. PYK_4dpi vs PYK
# A. PC
# pyk_4dpi vs pyk within the PC cluster
# first, looking at the pyk-4dpi induced genes
pc_dge_genes_pyk4dpi_pyk <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                    ident.1 = "pyk_4dpi", 
                                    ident.2 = "pyk_naive", 
                                    subset.ident = "PC", min.pct = 0.25,
                                    only.pos = F)
head(pc_dge_genes_pyk4dpi_pyk, n = 20)
dim(pc_dge_genes_pyk4dpi_pyk)



write.csv(pc_dge_genes_pyk4dpi_pyk, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/pc_dge_genes_pyk4dpi_pyk_081522_pm.csv')

EnhancedVolcano(pc_dge_genes_pyk4dpi_pyk,
                lab = rownames(pc_dge_genes_pyk4dpi_pyk),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c("H2-Q7", "Gbp4", "Stat1", "Ifitm3", "Bst2", "Gbp6",
                              "Gbp7", "Ada", "B2m", "H2-D1", "H2-K1", "H2-T23", 
                              "H2-Q6", "Ccnd", "Ero1b", "Mrpl57", "Rpl41", "Ly6c2",
                              "Nars", "Spcs3", "Gadd45gip1", "Lars2", "Jchain",
                              "B2m", "Enpp1", "Slamf9", "Zbp1", "Rbm3", "Ptprcap"),
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5,
                xlim = c(-3, 3),
                ylim = c(-0.1, 10),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# 02. PYK_4dpi vs PYK
# B. GC
# pyk_4dpi vs pyk within the PC cluster
# first, looking at the pyk-4dpi induced genes
gc_dge_genes_pyk4dpi_pyk <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                        ident.1 = "pyk_4dpi", 
                                        ident.2 = "pyk_naive", 
                                        subset.ident = "GC", min.pct = 0.25,
                                        only.pos = F)
head(gc_dge_genes_pyk4dpi_pyk, n = 20)
dim(gc_dge_genes_pyk4dpi_pyk)

write.csv(gc_dge_genes_pyk4dpi_pyk, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_pyk_081522_pm.csv')

EnhancedVolcano(gc_dge_genes_pyk4dpi_pyk,
                lab = rownames(gc_dge_genes_pyk4dpi_pyk),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c("Pdia4", "Calr", "Cct6a", "Fkbp2", "Pdia3", "Zbp1",
                              "Gbp2", "H2-Q7", "Irgm1", "Irf1", "Gbp4", "Stat1", 
                              "Bst2", "Gbp6", "Gbp7", "Ly6a",
                              "Xbp1", "H2-K1", "H2-Q7", "Irf1", "Ighg2b", 
                              "Iglc1", "Ighg2c", "Cd79a", "Cd79b", "Siglecg", "Rbm3", "Ptprcap"),
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 25),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# Draw a heatmap of B cells subset
#top50 <- findallmarkers %>% group_by(cluster) %>% top_n(n = 50)
#DoHeatmap(all5_dr32, features = top50$gene, group.bar = TRUE) + scale_fill_gradientn(colors = c("blue", "white", "red"))

# drawing venn diagram
venn_pky_adj_gc <- draw.pairwise.venn(area1 = 32, 
                                      area2 = 4, 
                                      cross.area = 920,
                                      category = c("Adjuvant" , "TbPYK"),
                                      fill = c("light blue", "pink"), 
                                      alpha = rep(0.5, 2), scaled=TRUE,
                                      cat.pos = c(0, 0),
                                      main = "TbPYK vs adj (GC)",
                                      fontfamily = rep("serif", 3))
grid.draw(venn_pky_adj_gc)


# drawing venn diagram
venn_pky_adj_pc <- draw.pairwise.venn(area1 = 5, 
                                      area2 = 0, 
                                      cross.area = 167,
                                      category = c("Adjuvant" , "TbPYK"),
                                      fill = c("light blue", "pink"), 
                                      main = "TbPYK vs adj (PC)",
                                      fontfamily = rep("serif", 3))
grid.draw(venn_pky_adj_pc)

# 03. PYK_4dpi vs ADJ_4dpi (what is actually reducing the parasitemia?)
# A. GC
# PYK_4dpi vs ADJ_4dpi within the GC cluster
# first, looking at the pyk-4dpi induced genes
gc_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                        ident.1 = "pyk_4dpi", 
                                        ident.2 = "adj_4dpi", 
                                        subset.ident = "GC", min.pct = 0.25,
                                        only.pos = F)
head(gc_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(gc_dge_genes_pyk4dpi_adj_4dpi)

write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(gc_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(gc_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 16),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# 03. PYK_4dpi vs ADJ_4dpi (what is actually reducing the parasitemia?)
# B. PC
# PYK_4dpi vs ADJ_4dpi within the PC cluster
# first, looking at the pyk-4dpi induced genes
pc_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                             ident.1 = "pyk_4dpi", 
                                             ident.2 = "adj_4dpi", 
                                             subset.ident = "PC", min.pct = 0.25,
                                             only.pos = F)
head(pc_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(pc_dge_genes_pyk4dpi_adj_4dpi)

write.csv(pc_dge_genes_pyk4dpi_adj_4dpi, 
          file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/pc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(pc_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(pc_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                selectLab = c("Lars2", "Ccnd2", "Tmsb4x", "Rpl22l1", "Actg1",
                              "H2-T23", "Actg1", "Rbm3", "Dph3", "Ptprcap",
                              "Snrpf"),
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 6),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# 03. PYK_4dpi vs ADJ_4dpi (what is actually reducing the parasitemia?)
# C. MzB
# PYK_4dpi vs ADJ_4dpi within the MzB cluster
# first, looking at the pyk-4dpi induced genes
mzb_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                             ident.1 = "pyk_4dpi", 
                                             ident.2 = "adj_4dpi", 
                                             subset.ident = "MzB", min.pct = 0.25,
                                             only.pos = F)
head(mzb_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(mzb_dge_genes_pyk4dpi_adj_4dpi)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
          #file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(mzb_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(mzb_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 100),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# D. FoB
# PYK_4dpi vs ADJ_4dpi within the FoB cluster
# first, looking at the pyk-4dpi induced genes
fob_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "adj_4dpi", 
                                              subset.ident = "FoB", min.pct = 0.25,
                                              only.pos = F)
head(fob_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(fob_dge_genes_pyk4dpi_adj_4dpi)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(fob_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(fob_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 100),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# E. B1 B
# PYK_4dpi vs ADJ_4dpi within the B1 B cluster
# first, looking at the pyk-4dpi induced genes
b1b_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "adj_4dpi", 
                                              subset.ident = "B1 B", min.pct = 0.25,
                                              only.pos = F)
head(b1b_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(b1b_dge_genes_pyk4dpi_adj_4dpi)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(b1b_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(b1b_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-1.5, 1.5),
                ylim = c(-0.1, 5),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# F. T1 B
# PYK_4dpi vs ADJ_4dpi within the T1 B cluster
# first, looking at the pyk-4dpi induced genes
t1b_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "adj_4dpi", 
                                              subset.ident = "T1 B", min.pct = 0.25,
                                              only.pos = F)
head(t1b_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(t1b_dge_genes_pyk4dpi_adj_4dpi)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(t1b_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(t1b_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 16),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# G. T2/T3 B
# PYK_4dpi vs ADJ_4dpi within the T2/T3 B cluster
# first, looking at the pyk-4dpi induced genes
t2t3b_dge_genes_pyk4dpi_adj_4dpi <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "adj_4dpi", 
                                              subset.ident = "T2/T3 B", min.pct = 0.25,
                                              only.pos = F)
head(t2t3b_dge_genes_pyk4dpi_adj_4dpi, n = 20)
dim(t2t3b_dge_genes_pyk4dpi_adj_4dpi)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(t2t3b_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(t2t3b_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2, 2),
                ylim = c(-0.1, 10),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")



# 03. PYK_4dpi vs PYK (what is actually reducing the parasitemia?)
# C. MzB
# PYK_4dpi vs pyk within the MzB cluster
# first, looking at the pyk-4dpi induced genes
mzb_dge_genes_pyk4dpi_pyk <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "pyk_naive", 
                                              subset.ident = "MzB", min.pct = 0.25,
                                              only.pos = F)
head(mzb_dge_genes_pyk4dpi_pyk, n = 20)
dim(mzb_dge_genes_pyk4dpi_pyk)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(mzb_dge_genes_pyk4dpi_pyk,
                lab = rownames(mzb_dge_genes_pyk4dpi_pyk),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-1, 1),
                ylim = c(-0.1, 16),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# D. FoB
# PYK_4dpi vs pyk within within the FoB cluster
# first, looking at the pyk-4dpi induced genes
fob_dge_genes_pyk4dpi_pyk <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                              ident.1 = "pyk_4dpi", 
                                              ident.2 = "pyk_naive", 
                                              subset.ident = "FoB", min.pct = 0.25,
                                              only.pos = F)
head(fob_dge_genes_pyk4dpi_pyk, n = 20)
dim(fob_dge_genes_pyk4dpi_pyk)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(fob_dge_genes_pyk4dpi_adj_4dpi,
                lab = rownames(fob_dge_genes_pyk4dpi_adj_4dpi),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-3, 3),
                ylim = c(-0.1, 16),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")



# S1. FoB - infection-induced
# ADJ_4dpi vs naive within within the FoB cluster
# first, looking at the pyk-4dpi induced genes
fob_dge_genes_adj4dpi_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                         ident.1 = "adj_4dpi", 
                                         ident.2 = "naive", 
                                         subset.ident = "FoB", min.pct = 0.25,
                                         only.pos = F)
head(fob_dge_genes_adj4dpi_naive, n = 20)
dim(fob_dge_genes_adj4dpi_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(fob_dge_genes_adj4dpi_naive,
                lab = rownames(fob_dge_genes_adj4dpi_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-12,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-5, 5),
                ylim = c(-0.1, 120),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")

# S2. MzB - infection-induced
# ADJ_4dpi vs naive within within the MzB cluster
# first, looking at the pyk-4dpi induced genes
mzb_dge_genes_adj4dpi_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                           ident.1 = "adj_4dpi", 
                                           ident.2 = "naive", 
                                           subset.ident = "MzB", min.pct = 0.25,
                                           only.pos = F)
head(mzb_dge_genes_adj4dpi_naive, n = 20)
dim(mzb_dge_genes_adj4dpi_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(mzb_dge_genes_adj4dpi_naive,
                lab = rownames(mzb_dge_genes_adj4dpi_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-12,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2.5, 2.5),
                ylim = c(-0.1, 120),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# S3. MzB - naive vs adj
# ADJ_4dpi vs naive within within the MzB cluster
# first, looking at the pyk-4dpi induced genes
mzb_dge_genes_adj_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                           ident.1 = "adj_naive", 
                                           ident.2 = "naive", 
                                           subset.ident = "MzB", min.pct = 0.25,
                                           only.pos = F)
head(mzb_dge_genes_adj_naive, n = 20)
dim(mzb_dge_genes_adj_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(mzb_dge_genes_adj_naive,
                lab = rownames(mzb_dge_genes_adj_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2.5, 2.5),
                ylim = c(-0.1, 120),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# S4. FoB - naive vs adj
# ADJ_4dpi vs naive within within the MzB cluster
# first, looking at the pyk-4dpi induced genes
fob_dge_genes_adj_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                       ident.1 = "adj_naive", 
                                       ident.2 = "naive", 
                                       subset.ident = "FoB", min.pct = 0.25,
                                       only.pos = F)
head(fob_dge_genes_adj_naive, n = 20)
dim(fob_dge_genes_adj_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(fob_dge_genes_adj_naive,
                lab = rownames(fob_dge_genes_adj_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2.5, 2.5),
                ylim = c(-0.1, 120),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# S5. GC - naive vs adj
# ADJ_4dpi vs naive within within the MzB cluster
# first, looking at the pyk-4dpi induced genes
gc_dge_genes_adj_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                       ident.1 = "adj_naive", 
                                       ident.2 = "naive", 
                                       subset.ident = "GC", min.pct = 0.25,
                                       only.pos = F)
head(gc_dge_genes_adj_naive, n = 20)
dim(gc_dge_genes_adj_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(gc_dge_genes_adj_naive,
                lab = rownames(gc_dge_genes_adj_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2.5, 2.5),
                ylim = c(-0.1, 10),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


# S6. PC - naive vs adj
# ADJ_4dpi vs naive within within the MzB cluster
# first, looking at the pyk-4dpi induced genes
pc_dge_genes_adj_naive <- FindMarkers(bcell_subset_final, group.by="orig.ident", 
                                      ident.1 = "adj_naive", 
                                      ident.2 = "naive", 
                                      subset.ident = "PC", min.pct = 0.25,
                                      only.pos = F)
head(pc_dge_genes_adj_naive, n = 20)
dim(pc_dge_genes_adj_naive)

#write.csv(gc_dge_genes_pyk4dpi_adj_4dpi$avg_log2FC>0.5, 
#file='/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_dge_genes_pyk4dpi_adj_4dpi_081522_pm.csv')

EnhancedVolcano(pc_dge_genes_adj_naive,
                lab = rownames(pc_dge_genes_adj_naive),
                x = 'avg_log2FC',
                y = 'p_val',
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                pCutoff = 10e-4,
                FCcutoff = 0.5, # this is a 2-fold change
                xlim = c(-2.5, 2.5),
                ylim = c(-0.1, 10),
                pointSize = 2.0,
                labSize = 3.0,
                title = NULL,
                subtitle = NULL,
                caption = NULL,
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 1,
                legendPosition = 'right',
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black') + theme(legend.position="none")


#######################################################################
# 03_new_Gene ontology analysis 081522_pm
#######################################################################

# 01. PYK vs ADJ
# A. PC
# pyk_naive vs adj_naive within the PC cluster
# first, looking at the pyk-vaccinated induced genes


# TRYING GSEA SCRIPTS
# Not working
msigdbr_species()
m_t2g <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)

em.pc_pyk_adj <- GSEA(df_pc_dge_genes_pyk_adj_gene_list, 
                      TERM2GENE = m_t2g,
                      verbose = FALSE,
                      exponent = 1,
                      minGSSize = 10,
                      maxGSSize = 500,
                      eps = 1e-10,
                      pAdjustMethod = "BH")

em.pc_pyk_adj[1:5,1:5]
doplot(em.pc_pyk_adj,showCategory=12,split=".sign")+facet_grid(~.sign)
gseaplot2(em.pc_pyk_adj, geneSetID = "HALLMARK_TNFA_SIGNALING_VIA_NFKB", title = em.cluster8$ID["HALLMARK_TNFA_SIGNALING_VIA_NFKB"],color="red",pvalue_table=T)                                                     
gseaplot2(em.pc_pyk_adj, geneSetID = c(1,2,3,5), title ="Activated Pathway",pvalue_table=T)

####

#new_pc_dge_genes_pyk_adj <- read.csv("/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/gc_pyk_dge_manual.csv")
#head(new_pc_dge_genes_pyk_adj)

df_pc_pyk_adj <- pc_dge_genes_pyk_adj
df_pc_pyk_adj

# filter on min log2fold change (log2FoldChange > 2)
gene_list_pc_pyk_adj <- df_pc_pyk_adj$avg_log2FC>0.5  
gene_list_pc_pyk_adj <- df_pc_pyk_adj$p_val<0.001

head(gene_list_pc_pyk_adj, n=10)

# name the vector
names(gene_list_pc_pyk_adj) <- as.character(row.names(df_pc_pyk_adj))
type(gene_list_pc_pyk_adj)
head(gene_list_pc_pyk_adj, n=10)

# omit any NA values 
# sort the list in decreasing order (required for clusterProfiler)
gene_list_pc_pyk_adj = sort(gene_list_pc_pyk_adj, decreasing = TRUE)
gene_list_pc_pyk_adj

df_gene_list_pc_pyk_adj <- names(gene_list_pc_pyk_adj)
df_gene_list_pc_pyk_adj


### TRIAL AGAIN ###

# we want the log2 fold change 
original_gene_list <- pc_pyk_adj$log2FoldChange

# name the vector
names(original_gene_list) <- pc_pyk_adj$gene

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.001)
sig_genes_df = subset(pc_pyk_adj, p_val < 0.001)

# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$avg_log2FC

# Name the vector
names(genes) <- sig_genes_df$gene

# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]


#-----------Over representation GO Analysis---------------

# ALL
ego_pc_all <- enrichGO(gene          = df_gene_list_pc_pyk_adj,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = "SYMBOL",
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

# BP (Biological Process)
ego_pc_bp <- enrichGO(gene          = df_gene_list_pc_pyk_adj,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "SYMBOL",
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.001,
                      qvalueCutoff  = 0.05)

head(ego_pc_bp)

# MF (Molecular Function)
ego_pc_mf <- enrichGO(gene          = df_gene_list_pc_pyk_adj,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "SYMBOL",
                      ont           = "MF",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

# CC (Cellular Component)
ego_pc_cc <- enrichGO(gene          = df_gene_list_pc_pyk_adj,
                      OrgDb         = org.Mm.eg.db,
                      keyType       = "SYMBOL",
                      ont           = "CC",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05)

barplot(ego_pc_all,showCategory=20, color='pvalue')
barplot(ego_pc_bp,showCategory=20, color='pvalue', title = "PYK vs adj (PC)")
barplot(ego_pc_mf,showCategory=20, color='pvalue', title = "PYK vs adj (PC)")
barplot(ego_pc_cc,showCategory=20, color='pvalue')

dotplot(ego_pc_all)
dotplot(ego_pc_bp, title = "PYK vs adj (PC)")
dotplot(ego_pc_mf, title = "PYK vs adj (PC)")
dotplot(ego_pc_cc)

cnetplot(ego_pc_all, categorySize="pvalue", foldChange = gene_list_pc, circular=FALSE)
cnetplot(ego_pc_all, foldChange=gene_list_pc, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_pc_bp, categorySize="pvalue", foldChange = gene_list_pc, circular=FALSE)
cnetplot(ego_pc_bp, foldChange=gene_list_pc, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_pc_mf, categorySize="pvalue", foldChange = gene_list_pc, circular=FALSE)
cnetplot(ego_pc_mf, foldChange=gene_list_pc, circular = TRUE, colorEdge = TRUE)


heatplot(ego_pc_all, showCategory = 10, foldChange = gene_list_pc)
heatplot(ego_pc_bp, showCategory = 20, foldChange = gene_list_pc)
heatplot(ego_pc_mf, showCategory = 20, foldChange = gene_list_pc)


ego_pc_all_pair <- pairwise_termsim(ego_pc_all)
emapplot(ego_pc_all_pair, color='pvalue')


# B. GC
# pyk_naive vs adj_naive within the GC cluster
# first, looking at the pyk-vaccinated induced genes

df_gc_pyk_adj <- gc_dge_genes_pyk_adj
df_gc_pyk_adj

# filter on min log2fold change (log2FoldChange > 2)
gene_list_gc_pyk_adj <- df_gc_pyk_adj$avg_log2FC>0.5  
gene_list_gc_pyk_adj <- df_gc_pyk_adj$p_val<0.001

# name the vector
names(gene_list_gc_pyk_adj) <- as.character(row.names(df_gc_pyk_adj))
type(gene_list_gc_pyk_adj)

# omit any NA values 
# sort the list in decreasing order (required for clusterProfiler)
gene_list_gc_pyk_adj = sort(gene_list_gc_pyk_adj, decreasing = TRUE)
gene_list_gc_pyk_adj

df_gene_list_gc_pyk_adj <- names(gene_list_gc_pyk_adj)
df_gene_list_gc_pyk_adj


#-----------Over representation GO Analysis---------------

# ALL
ego_gc_pyk_adj_all <- enrichGO(gene          = df_gene_list_gc_pyk_adj,
                               OrgDb         = org.Mm.eg.db,
                               keyType       = "SYMBOL",
                               ont           = "ALL",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05)

# BP (Biological Process)
ego_gc_pyk_adj_bp <- enrichGO(gene          = df_gene_list_gc_pyk_adj,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.001,
                              qvalueCutoff  = 0.05)

# MF (Molecular Function)
ego_gc_pyk_adj_mf <- enrichGO(gene          = df_gene_list_gc_pyk_adj,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "MF",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.001,
                              qvalueCutoff  = 0.05)

# CC (Cellular Component)
ego_gc_pyk_adj_cc <- enrichGO(gene          = df_gene_list_gc_pyk_adj,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "CC",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.001,
                              qvalueCutoff  = 0.05)

barplot(ego_gc_pyk_adj_all,showCategory=20, color='pvalue')
barplot(ego_gc_pyk_adj_bp,showCategory=20, color='pvalue', title = "PYK vs adj (GC)")
barplot(ego_gc_pyk_adj_mf,showCategory=20, color='pvalue')
barplot(ego_gc_pyk_adj_cc,showCategory=20, color='pvalue')

dotplot(ego_gc_pyk_adj_all)
dotplot(ego_gc_pyk_adj_bp, title = "PYK vs adj (GC)")
dotplot(ego_gc_pyk_adj_mf, title = "PYK vs adj (GC)")
dotplot(ego_gc_pyk_adj_cc)

cnetplot(ego_gc_pyk_adj_all, categorySize="pvalue", foldChange = gene_list_gc_pyk_adj, circular=FALSE)
cnetplot(ego_gc_pyk_adj_all, foldChange=gene_list_gc_pyk_adj, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_gc_pyk_adj_bp, categorySize="pvalue", foldChange = gene_list_gc_pyk_adj, circular=FALSE)
cnetplot(ego_gc_pyk_adj_bp, foldChange=gene_list_gc_pyk_adj, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_gc_pyk_adj_mf, categorySize="pvalue", foldChange = gene_list_gc_pyk_adj, circular=FALSE)
cnetplot(ego_gc_pyk_adj_mf, foldChange=gene_list_gc_pyk_adj, circular = TRUE, colorEdge = TRUE)

heatplot(ego_gc_pyk_adj_all, showCategory = 10, foldChange = gene_list_gc_pyk_adj)
heatplot(ego_gc_pyk_adj_bp, showCategory = 20, foldChange = gene_list_gc_pyk_adj)
heatplot(ego_gc_pyk_adj_mf, showCategory = 20, foldChange = gene_list_gc_pyk_adj)


ego_gc_pyk_adj_all_pair <- pairwise_termsim(ego_gc_pyk_adj_all)
emapplot(ego_gc_pyk_adj_all_pair, color='pvalue')


# 02. PYK_4dpi vs PYK
# A. PC
# pyk_4dpi vs pyk within the PC cluster
# first, looking at the pyk-4dpi induced genes

df_pc_pyk4dpi_pyk <- pc_dge_genes_pyk4dpi_pyk
df_pc_pyk4dpi_pyk

# filter on min log2fold change (log2FoldChange > 2)
gene_list_pc_pyk4dpi_pyk <- df_pc_pyk4dpi_pyk$avg_log2FC>0.5
gene_list_pc_pyk4dpi_pyk <- df_pc_pyk4dpi_pyk$p_val<0.001

# name the vector
names(gene_list_pc_pyk4dpi_pyk) <- as.character(row.names(df_pc_pyk4dpi_pyk))
type(gene_list_pc_pyk4dpi_pyk)

# omit any NA values 
# sort the list in decreasing order (required for clusterProfiler)
gene_list_pc_pyk4dpi_pyk = sort(gene_list_pc_pyk4dpi_pyk, decreasing = TRUE)
gene_list_pc_pyk4dpi_pyk

df_pc_pyk4dpi_pyk <- names(gene_list_pc_pyk4dpi_pyk)
df_pc_pyk4dpi_pyk


#-----------Over representation GO Analysis---------------

# ALL
ego_pc_pyk4dpi_pyk_all <- enrichGO(gene          = df_pc_pyk4dpi_pyk,
                                   OrgDb         = org.Mm.eg.db,
                                   keyType       = "SYMBOL",
                                   ont           = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05)

# BP (Biological Process)
ego_pc_pyk4dpi_pyk_bp <- enrichGO(gene          = df_pc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.001,
                                  qvalueCutoff  = 0.05)

head(ego_pc_pyk4dpi_pyk_bp)

# MF (Molecular Function)
ego_pc_pyk4dpi_pyk_mf <- enrichGO(gene          = df_pc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "MF",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05)

# CC (Cellular Component)
ego_pc_pyk4dpi_pyk_cc <- enrichGO(gene          = df_pc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "CC",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05)

barplot(ego_pc_pyk4dpi_pyk_all,showCategory=20, color='pvalue')
barplot(ego_pc_pyk4dpi_pyk_bp,showCategory=20, color='pvalue', title = "PYK vs adj (PC)")
barplot(ego_pc_pyk4dpi_pyk_mf,showCategory=20, color='pvalue')
barplot(ego_pc_pyk4dpi_pyk_cc,showCategory=20, color='pvalue')

dotplot(ego_pc_pyk4dpi_pyk_all)
dotplot(ego_pc_pyk4dpi_pyk_bp, title = "PYK+4dpi vs PYK (PC)")
dotplot(ego_pc_pyk4dpi_pyk_mf)
dotplot(ego_pc_pyk4dpi_pyk_cc)

cnetplot(ego_pc_pyk4dpi_pyk_all, categorySize="pvalue", foldChange = df_pc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_pc_pyk4dpi_pyk_all, foldChange=df_pc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_pc_pyk4dpi_pyk_bp, categorySize="pvalue", foldChange = df_pc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_pc_pyk4dpi_pyk_bp, foldChange=df_pc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_pc_pyk4dpi_pyk_mf, categorySize="pvalue", foldChange = df_pc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_pc_pyk4dpi_pyk_mf, foldChange=df_pc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

heatplot(ego_pc_pyk4dpi_pyk_all, showCategory = 10, foldChange = df_pc_pyk4dpi_pyk)
heatplot(ego_pc_pyk4dpi_pyk_bp, showCategory = 20, foldChange = df_pc_pyk4dpi_pyk)
heatplot(ego_pc_pyk4dpi_pyk_mf, showCategory = 20, foldChange = df_pc_pyk4dpi_pyk)


ego_pc_pyk4dpi_pyk_all_pair <- pairwise_termsim(ego_pc_pyk4dpi_pyk_all)
emapplot(ego_pc_pyk4dpi_pyk_all_pair, color='pvalue')



# B. GC
# pyk_4dpi vs pyk within the GC cluster
# first, looking at the pyk-4dpi induced genes

df_gc_pyk4dpi_pyk <- gc_dge_genes_pyk4dpi_pyk
df_gc_pyk4dpi_pyk

# filter on min log2fold change (log2FoldChange > 2)
gene_list_gc_pyk4dpi_pyk <- df_gc_pyk4dpi_pyk$avg_log2FC>0.5

# name the vector
names(gene_list_gc_pyk4dpi_pyk) <- as.character(row.names(df_gc_pyk4dpi_pyk))
type(gene_list_gc_pyk4dpi_pyk)

# omit any NA values 
# sort the list in decreasing order (required for clusterProfiler)
gene_list_gc_pyk4dpi_pyk = sort(gene_list_gc_pyk4dpi_pyk, decreasing = TRUE)
gene_list_gc_pyk4dpi_pyk

df_gc_pyk4dpi_pyk <- names(gene_list_gc_pyk4dpi_pyk)
df_gc_pyk4dpi_pyk


#-----------Over representation GO Analysis---------------

# ALL
ego_gc_pyk4dpi_pyk_all <- enrichGO(gene          = df_gc_pyk4dpi_pyk,
                                   OrgDb         = org.Mm.eg.db,
                                   keyType       = "SYMBOL",
                                   ont           = "ALL",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.001,
                                   qvalueCutoff  = 0.05)

# BP (Biological Process)
ego_gc_pyk4dpi_pyk_bp <- enrichGO(gene          = df_gc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "BP",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.001,
                                  qvalueCutoff  = 0.05)

# MF (Molecular Function)
ego_gc_pyk4dpi_pyk_mf <- enrichGO(gene          = df_gc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "MF",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05)

# CC (Cellular Component)
ego_gc_pyk4dpi_pyk_cc <- enrichGO(gene          = df_gc_pyk4dpi_pyk,
                                  OrgDb         = org.Mm.eg.db,
                                  keyType       = "SYMBOL",
                                  ont           = "CC",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05)

barplot(ego_gc_pyk4dpi_pyk_all,showCategory=20, color='pvalue')
barplot(ego_gc_pyk4dpi_pyk_bp,showCategory=20, color='pvalue', title = "PYK+4dpi vs PYK (GC)")
barplot(ego_gc_pyk4dpi_pyk_mf,showCategory=20, color='pvalue')
barplot(ego_gc_pyk4dpi_pyk_cc,showCategory=20, color='pvalue')

dotplot(ego_gc_pyk4dpi_pyk_all)
dotplot(ego_gc_pyk4dpi_pyk_bp, title = "PYK+4dpi vs PYK (GC)")
dotplot(ego_gc_pyk4dpi_pyk_mf, title = "PYK+4dpi vs PYK (GC)")
dotplot(ego_gc_pyk4dpi_pyk_cc)

cnetplot(ego_gc_pyk4dpi_pyk_all, categorySize="pvalue", foldChange = df_gc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_gc_pyk4dpi_pyk_all, foldChange=df_gc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_gc_pyk4dpi_pyk_bp, categorySize="pvalue", foldChange = df_gc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_gc_pyk4dpi_pyk_bp, foldChange=df_gc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

cnetplot(ego_gc_pyk4dpi_pyk_mf, categorySize="pvalue", foldChange = df_gc_pyk4dpi_pyk, circular=FALSE)
cnetplot(ego_gc_pyk4dpi_pyk_mf, foldChange=df_gc_pyk4dpi_pyk, circular = TRUE, colorEdge = TRUE)

heatplot(ego_gc_pyk4dpi_pyk_all, showCategory = 10, foldChange = df_gc_pyk4dpi_pyk)
heatplot(ego_gc_pyk4dpi_pyk_bp, showCategory = 20, foldChange = df_gc_pyk4dpi_pyk)
heatplot(ego_gc_pyk4dpi_pyk_mf, showCategory = 20, foldChange = gene_list_gc_pyk4dpi_pyk)


ego_gc_pyk4dpi_pyk_all_pair <- pairwise_termsim(ego_gc_pyk4dpi_pyk_all)
emapplot(ego_gc_pyk4dpi_pyk_all_pair, color='pvalue')



'''
#######################################################################
# 02_new DGE
# As a default, Seurat uses the non-parametric Wilcoxon rank sum test. 
#######################################################################

# read previously pre-processed, dimensional analysis, and clustered file
all5_bcellsubset <- readRDS(file="/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/02_all5_res75_dr32_bcellsubset_080822_pm.rds")

# Naming clusters based on DEG
new_cluster_names <- c("FoB", "T1 B", "T2/T3 B", "MzB", "B1 B/Breg",
                       "GC", "Pro B", "Usp18+, Ifit+", "MBC", "9_macro", "10", 
                       "PC", "12_nonB", "13_nonB")
names(new_cluster_names) <- levels(all5_bcellsubset)
all5_bcellsubset <- RenameIdents(all5_bcellsubset, new_cluster_names)

# remove non-B cells or multiplets
bcellsubset_removal <- subset((all5_bcellsubset), 
                              idents = c("MBC", "9_macro", "10", "12_nonB", "13_nonB", "Pro B"), 
                              invert = T)
DimPlot(bcellsubset_removal, label = T, pt.size = 0.1)

#DimPlot(all5_bcellsubset, reduction = "umap",  label = T, pt.size = 0.1)
#DimPlot(all5_bcellsubset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.1)
#DimPlot(all5_bcellsubset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.1)


# How many cells are in each cluster
table(Idents(bcellsubset_removal))
all5_annotation
head(bcellsubset_removal@meta.data)

### Cell number analysis ###
# How many cells are in each cluster
B_cell_counts <- table(Idents(bcellsubset_removal), bcellsubset_removal$group)
head(B_cell_counts, n =15)

# filter out non-B cells

# cluster 11
dc <- c("Itgam", "Itgax", "Il12a", "Ifng", "Trbc1", "Ccr7", "Fcer1a", "Cst3")
FeaturePlot(all5_bcellsubset, features = dc, label = F)
DotPlot(all5_bcellsubset, features = dc, cols = c("blue", "red"))

# cluster 11
tcell <- c("Trac", "Trbc1", "Cd3e", "Cd4", "Cd8a", "Foxp3") 
FeaturePlot(all5_bcellsubset, features = tcell, label = F)

# Macrophages
# cluster 9
FeaturePlot(all5_bcellsubset, features = c("C1qa", "Ms4a7"), label = F) 
# NK cells 
# cluster 11
FeaturePlot(all5_bcellsubset, features = c("Ncr1", "Klrb1c", "Nkg7"), label = F)
# Neutrophils/Monocytes
FeaturePlot(all5_bcellsubset, features = c("Ly6c2", "Ear2", "Fcgr3a", "Ms4a7"), label = F)



# identifying each cluster with well-known canonical markers

c14_genes <- c("Cd19", "Ms4a1", "Cd40", "Ighd", "Prdm1", "Tnfrsf13b",
               "Fcer2a", "Pou2af1", "Aicda", "Igha", "Ighg1",
               "Ighg2b", "Ighg3", "Sdc1", "Mki67", "Ly6d", "Ly6c2",
               "Cd3d", "Icos", "Cxcr3", "Ifng", "Ccr6", "Il22", "Pdcd1", 
               "Pdcd1lg2", "Cd44",
               "Cd4", "Cd8a", "Foxp3", "Cd3e", "Trdc", "Il7r", "Ncr1", "Kit", 
               "Ms4a2", "Lyve1", "Itgax", "Cxcr6", "Ccr7", 
               "Xcr1", "Clec9a", "Clec10a", "Irf4", "Cd14",  
               "Bcl6", "Cd1d1", "Cd9", "Pax5", "Ptprc", "Xbp1", "Jchain",
               "Cd38", "Hhex", "Zeb2", "Klf2", "Zbtb32", "Bach2",
               "Sell", "Id2", "Klrd1", "S1pr1", "Il10", "Il10ra",
               "C1qa", "Ear2", "Ighm", "Iglc1", "Iglc2",
               "Prss2", "Zg16", "Ctrl", "Reg1", "Reg2", "Rpl13",
               "Usp18", "Ifit2", "Isg20", "Cr2")

DotPlot(object = all5_bcellsubset, features = c14_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# prepro
prepro <- c("Ebf1", "Ybx3", "Cd93", "Il2ra", "Spn", "Il7r", "Flt3", "Kit", "Bok", "Cd300a")
FeaturePlot(all5_bcellsubset, features = prepro, label = T)

# cluster 13 and cluster 12 are non-B cells
# cluster 13 is multiplets showing Clec9a, Clec10a, Ly6d, Ly6c2, Cd4, Cd8a, Foxp3
# cluster 12 is T cell-stimulatory or DC showing Cd9
# also both cluster are very low in cell number

# cluster 10 is not a B cell named "Reg+, Rpl13+"
# possible a component of ribosome
# find all markers of cluster "Reg+, Rpl13+" vs rest
cluster10_markers <- FindMarkers(all5_bcellsubset, ident.1 = "Reg+, Rpl13+", only.pos = T, min.pct = 0.25)
head(cluster10_markers, n = 20)

# cluster 7 named "Usp18+, Ifit+"
# find all markers of cluster "Usp18+, Ifit+" vs rest
cluster7_markers <- FindMarkers(all5_bcellsubset, ident.1 = "Usp18+, Ifit+", only.pos = T, min.pct = 0.25)
head(cluster7_markers, n = 30)

# cluster 2 named "T2/T3 B"
# find all markers of cluster "T2/T3 B" vs rest
# T2 B cells = Cr2 (Cd21) absence & Fcer2a (Cd23), Cd24a, IgM and IgD high expression
# T3 B cells = IgD and Fcer2a (CD23) high, IgM and CD24a low, Cr2 (Cd21) absence
cluste2_markers <- FindMarkers(all5_bcellsubset, ident.1 = "T2/T3 B", only.pos = T, min.pct = 0.25)
head(cluster2_markers, n = 20)

c2_genes <- c("Ly6d", "Ms4a1", "Ighd", "Ighm", "Fcer2a", "Cd19", "Cr2", "Cd93",
              "Cd24a")
DotPlot(object = all5_bcellsubset, features = c2_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# cluster 9 is named "Fcer1g+, C1q+"
# find all markers of cluster "T3 B" vs rest
cluster9_markers <- FindMarkers(all5_bcellsubset, ident.1 = "Fcer1g+, C1q+", only.pos = T, min.pct = 0.25)
head(cluster9_markers, n = 20)

# check Memory B cells
# this could be T cells with expression of Cd3d, Cd8a
# cluster 9 or cluster 13?
# Id2, Ccr7, S1pr1, Klf2, S100a4, 
# Cd73 (Pdcd1), Cd273 (Pdcd1lg2), Itgax, Clec9a, Clec10a, Klrd1


mbc_genes <- c("Cd19", "Id2", "Hhex", "Zeb2",
               "Cd44", "Cxcl12",
               "Klrd1", "Pdcd1", "Pdcd1lg2", "Itgax", "Clec9a", "Clec10a",
               "Ccr7", "Cd40" ,"Klf2", "Zbtb32",
               "Ighd", "Tnfrsf13b", "Cxcr5", "Fcer2a",
               "Tagln2", "S100a4", "Itgb1",
               "Zap70", "Cd28", "Arl4c", "Cd247",
               "Cr2", "Cd52", "Vim", "Parm1", "Cd55", "Foxp1", "Cxcr4", # Riedel_nature_2020
               "Cxcr3", "Ptpn22", "Cd9",
               "Ighg1", "Ighg2b", "Ighg2c", "Ighg3",
               "Cd8a", "Foxp3", "Cd3e", "Cd4", # T cells?
               "Sdc1", "Mki67", "Ly6d", "Ly6c2", # NK?
               "S1pr1", "Il10", "Il10ra")

DotPlot(object = all5_bcellsubset, features = new_mbc_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")


# new mbc dot plot
new_mbc_genes <- c("Id2", "S100a4", "Pdcd1", "Pdcd1lg2", "Parm1", # Riedel_nature_2020
                   "Cd28", "Cd40", "Cd44", "Cd247", "Zap70", "Arl4c",
                   "Itgam", "Itgax", "Tbx21") #atMBC
DotPlot(object = bcell_subset_final, features = new_mbc_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

matureb <- c("Ms4a1", "Ly6a", "Fcer2a", "Bank1", "Bcl2", "Ltb", "Fcrl1", "Ctsh", "H2-Dm-a")
FeaturePlot(all5_bcellsubset, features = matureb, label = T)
DotPlot(all5_bcellsubset, features = matureb, cols = c("blue", "red"))

# B reg
breg <- c("Fcrl5", "Zbtb20", "Ccdc28b", "Cd9", "Ptpn22", "Atf3", "Il10", "Apoe")
FeaturePlot(all5_bcellsubset, features = breg, label = T)
DotPlot(object = all5_bcellsubset, features = breg) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

DimPlot(all5_bcellsubset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.1)
'''

### Cell number analysis ###
# How many cells are in each cluster
B_cell_counts <- table(Idents(bcellsubset_removal), bcellsubset_removal$group)
head(B_cell_counts, n =15)
