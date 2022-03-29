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
# COMMENTS
#############################################
# parameters for pre-processing are flexible.

# toc: Table of counts (filtered). Just those columns of \code{tod} that contain cells.
# tod: Table of droplets (raw).  A matrix with columns being each droplet and rows each gene.

# REMOVE BACKGROUND IN adj
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

# REMOVE BACKGROUND IN pyk
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

### Remove batch effects by Harmony
### iterates until similar cells are clustered using PCA for dimensionality reduction.

combine <- merge(naive, y = list(adj, pyk), add.cell.ids = c("naive", "adj", "pyk"), project = "PYK")
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
saveRDS(combine, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/harmony_integrated.rds")

#############################################
# 02 Single R Annotation
#############################################

ref <- ImmGenData()
ref

# Assigning general annotation
whole_annotation <- SingleR(test = SummarizedExperiment(combine), ref = ref, assay.type.test=1, 
                            labels = ref$label.main)

whole_annotation
head(combine@meta.data)
# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
combine$SingleR.pruned.calls <- whole_annotation$pruned.labels
combine$SingleR.calls <- whole_annotation$labels

# Run UMAP
combine <- RunUMAP(combine, dims = 1:30)
DimPlot(combine, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.5)
head(combine@meta.data)
combine

# Annotation Diagnostics
pred.grun <- SingleR(test=as.SingleCellExperiment(combine), ref=ref, assay.type.test=1, 
                     labels=ref$label.main, de.method="wilcox")
table(pred.grun$labels)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(combine) <- "SingleR.calls"

# How many cells are in each cluster
table(Idents(combine))

# How many cells are in each cluster
combine_counts <- table(Idents(combine))
head(combine_counts)

# to visualize the annotated table
combine_long <- as.data.frame(combine_counts)
head(combine_long)
colnames(combine_long) <- c("cell_type", "Cell_count")
head(combine_long)
rownames(combine_long) <- combine_long$x
combine_long

# Create ggplot2 plot scaled to 1.00
bar_counts <- ggplot(combine_long, aes(x = cell_type, y = Cell_count, fill = cell_type)) + 
  geom_bar(stat = "identity", color ="black") + scale_x_discrete(name ="Cell type") + 
  theme(legend.position="right") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Draw ggplot2 plot scaled to 1.00
bar_counts

# save annotated file
saveRDS(combine, file = "pyk_annotation_integrated_112121_PM.rds")

#############################################
# 03 B-cells Subsetting
#############################################

# Can I create a Seurat object of just the Plasma cells and B cells?
b_cells_general <- subset((combine))

# How many cells are in each cluster
table(Idents(b_cells_general))
b_cells_general
head(b_cells_general@meta.data)
# Visualize and just select B cells before subsetting.
DimPlot(b_cells_general, reduction = "umap",  label = T, pt.size = 1)

################**************################
# Now normalize B cells and cluster UMAP for subsetting.
b_cells_general <- NormalizeData(b_cells_general) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
ElbowPlot(b_cells_general, n=50)
b_cells_general <- b_cells_general %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
b_cells_general <- b_cells_general %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
# Visualize the newly normalized and clustered B cells
DimPlot(b_cells_general, reduction = "umap",  label = T, pt.size = 1)
################**************################

# Next, switch the identity class of all cells to reflect replicate ID
Idents(b_cells_general) <- "seurat_clusters"
results <- SingleR(test = as.SingleCellExperiment(b_cells_general), 
                   ref = ref, assay.type.test=1, 
                   labels = ref$label.main) # OR I could use label.main
results

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
b_cells_general$SingleR.pruned.calls_B <- results$pruned.labels
b_cells_general$SingleR.calls_B <- results$labels
head(b_cells_general@meta.data)

#Run UMAP
#b_cells_general <- RunUMAP(b_cells_general, dims = 1:30)
#DimPlot(b_cells_general, reduction = "umap", label = TRUE)
#DimPlot(b_cells_general, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
#DimPlot(b_cells_general, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
#DimPlot(b_cells_general, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(b_cells_general) <- "SingleR.calls_B"

# How many cells are in each cluster
table(Idents(b_cells_general))

# Can I create a Seurat object of just the Plasma cells and B cells?
b_cells_subset <- subset((b_cells_general), 
                         idents = c("B cells (B.CD19CONTROL)", "B cells (B.Fo)","B cells (B.FO)","B cells (B.FrE)", 
                                    "B cells (B.FRE)","B cells (B.FrF)","B cells (B.GC)","B cells (B.MEM)", 
                                    "B cells (B.MZ)","B cells (B.T1)","B cells (B.T2)","B cells (B.T3)", 
                                    "B cells (B1a)","B cells (B1A)","B cells (B1b)","B cells (preB.FrC)", 
                                    "B cells (preB.FrD)","B cells (preB.FRD)","B cells (proB.CLP)","B cells (proB.FrA)", 
                                    "B cells (proB.FrBC)"))
# Visualize the newly normalized and clustered B cells
DimPlot(b_cells_subset, reduction = "umap", label = T, pt.size = 1)

# How many cells are in each cluster
B_cell_counts <- table(Idents(b_cells_subset), b_cells_subset$group)
head(B_cell_counts)

B_cell_counts_long <- as.data.frame(B_cell_counts)    
colnames(B_cell_counts_long) <- c("B_cell", "Group", "Relative_abundance")
rownames(B_cell_counts_long) <- B_cell_counts_long$x
B_cell_counts_long

# Quantify different type of B cell genes
ggp <- ggplot(B_cell_counts_long, 
              aes(x = Group, 
                  y = Relative_abundance, 
                  fill = B_cell)) + 
  geom_bar(position = "fill", stat = "identity", color ="black") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_x_discrete(name ="Group") + theme(legend.position="right")
ggp  # Draw ggplot2 plot scaled to 1.00

# Save counts in CSV
#write.csv(B_cell_counts_long, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/Bcell_counts_112221_PM.txt")

################**************################
#Re-run UMAP on new created Bcell object for subsetting
b_cells_subset <- NormalizeData(b_cells_subset) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
ElbowPlot(b_cells_subset, n=50)
b_cells_subset <- b_cells_subset %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
b_cells_subset <- b_cells_subset %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(b_cells_subset, reduction = "umap",  label = T, pt.size = 1)
DimPlot(b_cells_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_subset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_subset, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)
################**************################

# Run UMAP
# b_cells_subset <- RunUMAP(b_cells_subset, dims = 1:30)
# DimPlot(b_cells_subset, reduction = "umap", label = TRUE)
# DimPlot(b_cells_subset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
# DimPlot(b_cells_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
# DimPlot(b_cells_subset, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

# DOT PLOT MARKER GENES
cd_genes <- c("Cd19", "Cd93", "Ly6d", "Ebf1","Ms4a1","Vpreb3","Itgam","Zbtb32","Zbtb20", 
              "Cr2","Ighd","Fcer2a","Aicda","Mki67","Ighm","Sdc1","Ezh2","Bcl6", "Cd1d1", 
              "Cd9","Pax5","Ptprc", "Xbp1", "Jchain")
DotPlot(object = b_cells_subset, features = cd_genes) + RotatedAxis() + 
  scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# B Cells Only
FeaturePlot(b_cells_subset, features = c("Cd79a", "Cd19", "Pax5", "Ebf1"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Cd79a", "Cd19", "Pax5", "Ebf1"), pt.size = 0.3, 
        group.by = "orig.ident")

# Transitional B cells
FeaturePlot(b_cells_subset, features = c("Ly6d", "Ebf1", "Ms4a1", "Vpreb3"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Ly6d", "Ebf1", "Ms4a1", "Vpreb3"), pt.size = 0.3, 
        group.by = "orig.ident") 

# Transitional T2 B cells
FeaturePlot(b_cells_subset, features = c("Ly6d", "Ms4a1", "Ighd", "Fcer2a"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Ly6d", "Ms4a1", "Ighd", "Fcer2a"), pt.size = 0.3, 
        group.by = "orig.ident") 

# Transitional T3 B cells
FeaturePlot(b_cells_subset, features = c("Cd19","Cr2","Ighd","Fcer2a"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Cd19", "Cr2", "Ighd", "Fcer2a"), pt.size = 0.3, 
        group.by = "orig.ident")

# B1A B cells
FeaturePlot(b_cells_subset, features = c("Zbtb20", "Cd19"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Zbtb20", "Cd19"), pt.size = 0.3, group.by = "orig.ident", 
        group.by = "orig.ident")

# B1B B cells
FeaturePlot(b_cells_subset, features = c("Zbtb32", "Itgam", "Cd19"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Zbtb32", "Itgam", "Cd19"), pt.size = 0.3, 
        group.by = "orig.ident")

# Atypical Memory B cells
FeaturePlot(b_cells_subset, features = c("Itgam","Zbtb32","Zbtb20", "Itgax", "Tbx21", "S100a6"), 
            pt.size = 0.3) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Itgam","Zbtb32","Zbtb20", "Itgax", "Tbx21", "S100a6"), 
        pt.size = 0.3, group.by = "orig.ident")

# Marginal Zone B cells
FeaturePlot(b_cells_subset, features = c("Ly6d", "Ebf1", "Cr2", "Ms4a1", "Cd19", "Cd1d1"), pt.size = 0.3) &
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Ly6d", "Ebf1", "Cr2", "Ms4a1", "Cd19", "Cd1d1"), pt.size = 0.3, 
        group.by = "orig.ident")

# Follicular B cells
FeaturePlot(b_cells_subset, features = c("Ighd", "Fcer2a", "Cd19", "Cr2", "Cd1d1"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Ighd", "Fcer2a", "Cd19", "Cr2", "Cd1d1"), pt.size = 0.3, 
        group.by = "orig.ident")

# Germinal center B cells
FeaturePlot(b_cells_subset, features = c("Aicda", "Mki67"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Aicda", "Mki67"), pt.size = 0.3, 
        group.by = "orig.ident")

# Plasma B cells
FeaturePlot(b_cells_subset, features = c("Sdc1", "Xbp1", "Ighm", "Jchain"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Sdc1", "Xbp1", "Ighm", "Jchain"), pt.size = 0.3, 
        group.by = "orig.ident")

### DIFFERENTIAL GENE EXPRESSION

b_cells_markers <- FindAllMarkers(b_cells_subset)
head(b_cells_markers)
write.csv(b_cells_markers, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/Bcell_subset_gex_022222_PM.txt")

EnhancedVolcano(b_cells_markers,
                lab = rownames(b_cells_markers),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 10e-4,
                FCcutoff = 1.5,
                xlim = c(-6, 6),
                pointSize = 2.0,
                labSize = 3.0,
                subtitle = NULL,
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-4',
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.5,)

levels(b_cells_subset)

# 1 trial
# Find differentially expressed features between 5 and 17 cluster
trial_1 <- FindMarkers(b_cells_subset, ident.1 = "5", ident.2 = "17")
# view results
head(trial_1)

# 2 trial for cluster 5
# Find differentially expressed features between cluster 5 and all other cells, only
# search for positive markers
cluster_5 <- FindMarkers(b_cells_subset, ident.1 = "5", ident.2 = NULL, only.pos = TRUE)
# view results
head(cluster_5)

# 2 trial for cluster 17
# Find differentially expressed features between cluster 17 and all other cells, only
# search for positive markers
cluster_17 <- FindMarkers(b_cells_subset, ident.1 = "17", ident.2 = NULL, only.pos = TRUE)
# view results
head(cluster_17)

# Find differentially expressed features between "B cells" and all other cells, 
# only search for positive markers
b_cells_specific_markers <- FindMarkers(b_cells_subset, ident.1 = "B cells", ident.2 = NULL, only.pos = TRUE)
# view results
head(b_cells_specific_markers)

# 1 trial -> I can use this information to uncover highly expressed genes in PC
# during infection stage at 4 dpi
# Find differentially expressed features between 5 and 17 cluster
# only search for positive markers
trial_1 <- FindMarkers(b_cells_subset, ident.1 = "5", ident.2 = "17", only.pos = TRUE)
# view results
head(trial_1)

# 3 trial for cluster 5
# Alternative DE test using MAST
# Find differentially expressed features between cluster 5 and all other cells,
# only search for positive markers
cluster_mast <- FindMarkers(b_cells_subset, ident.1 = "5", ident.2 = NULL, only.pos = TRUE, test.use = "MAST")
head(cluster_mast)

# Try this gene markers to annotate cell identities prior to ImmGen
FeaturePlot(combine, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# How do you know that these significant genes belong to which condition? 

## GENE REGULATORY NETWORKS
