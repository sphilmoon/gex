library(monocle3)
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

ref <- ImmGenData()

# Load annotated rds file
annotation_integrated <- readRDS(file = "/home/bmrc/Desktop/annotation_working.rds")

# 1st, filter: Subset T cells -----------------------------------------------
# Can I create a Seurat object of just the Plasma cells and B cells?
T_cells_1st <- subset((annotation_integrated), idents = c("T cells", "Tgd"))

# How many cells are in each cluster
table(Idents(T_cells_1st))

T_cells_1st
head(T_cells_1st@meta.data)
# Next, switch the identity class of all cells to reflect replicate ID
Idents(T_cells_1st) <- "seurat_clusters"


#Identify the B cell sub-clusters at a finer level using SingleR
results <- SingleR(test = as.SingleCellExperiment(T_cells_1st), ref = ref, assay.type.test=1,
                   labels = ref$label.fine)
results

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
T_cells_1st$SingleR.pruned.calls_B <- results$pruned.labels
T_cells_1st$SingleR.calls_B <- results$labels
head(T_cells_1st@meta.data)

# Now run umap to visualize the finer label from the 1st filter
T_cells_1st <- RunUMAP(T_cells_1st, dims = 1:10)
DimPlot(T_cells_1st, reduction = "umap", label = TRUE)
DimPlot(T_cells_1st, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cells_1st, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cells_1st, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)


# Next, switch the identity class of all cells to reflect replicate ID
Idents(T_cells_1st) <- "SingleR.calls_B"

# How many cells are in each cluster
table(Idents(T_cells_1st))

head(n=400, T_cells_1st)

# Type all the fine labels manually
# Is there other way to list these columns???
T_cell_subset <- subset((T_cells_1st), idents = c("T cells (T.4MEM49D+11A+.D30.LCMV)", "T cells (T.CD4TESTCJ)", 
                                                  "T cells (T.4EFF49D+11A+.D8.LCMV)", "T cells (T.4MEM44H62L)", "T cells (T.4MEM)",
                                                  "T cells (T.CD4.96H)", "T cells (T.CD4.48H)", 
                                                  "T cells (T.CD4TESTJS)", "T cells (T.4Mem)", "T cells (T.4MEM49D+11A+.D30.LCMV)",
                                                  "T cells (T.CD4.5H)", "T cells (T.4NVE44-49D-11A-)", "T cells (T.4SP24-)",
                                                  "T cells (T.CD4.24H)", "T cells (T.4Nve)",
                                                  "T cells (T.CD4.CTR)", "T cells (T.4.PLN)", "T cells (T.4FP3+25+)"))
                        
# How many cells are in each cluster
T_cell_counts <- table(Idents(T_cell_subset), T_cell_subset$group)
head(T_cell_counts)

data_long <- as.data.frame(T_cell_counts)    
colnames(data_long) <- c("T_cell", "Group", "Relative_abundance")

rownames(data_long) <- data_long$x

data_long

ggp <- ggplot(data_long,            # Create ggplot2 plot scaled to 1.00
              aes(x = Group,
                  y = Relative_abundance,
                  fill = T_cell)) +
  geom_bar(position = "fill", stat = "identity", color ="black") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  scale_x_discrete(name ="Group") +
  theme(legend.position="right")

ggp                                 # Draw ggplot2 plot scaled to 1.00
                        
#Run UMAP
T_cells_umap <- RunUMAP(T_cell_subset, dims = 1:10)
DimPlot(T_cell_subset, reduction = "umap", label = TRUE)
DimPlot(T_cell_subset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cell_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cell_subset, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)


# find markers for every cluster compared to all remaining cells, report only the positive ones
T_cell_markers <- FindAllMarkers(T_cell_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T_cell_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# list the gex CD4 T-cell markers
head(n=100, T_cell_markers)

# Use literature genes to define major CD4 gene. Top 15.
# drawing a dot plot
cd_genes <- c("Cd4", "Cd28", "Ifng", "Tnf", "Foxp3", "Ahr", "Ccr7", "Icos", "Il21", "Lag3", "Lrrc32", "Ccr5", "Cd44")
DotPlot(object = T_cell_subset, features = cd_genes) + RotatedAxis() + scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# T.4.PLN T cells genes
VlnPlot(T_cell_subset, features = c("Cd4", "Cd28", "Ifng", "Tnf", "Foxp3", "Ahr", "Ccr7", "Icos", "Il21", "Lag3", "Lrrc32", "Ccr5"), pt.size = 0.3, group.by = "orig.ident")  
FeaturePlot(T_cell_subset, features = c("Cd4", "Cd28", "Ifng", "Tnf", "Foxp3", "Ahr", "Ccr7", "Icos", "Il21", "Lag3", "Lrrc32", "Ccr5"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

# look at the most gex in my data, then compare to the literature
# Next, plot the most expressed genes in the CD4 subset.
# Naive
VlnPlot(T_cell_subset, features = c("Cd4", "Cd28", "Ccr7"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Cd4", "Cd28", "Ccr7"), pt.size = 0.3)

# Th17
VlnPlot(T_cell_subset, features = c("Il21", "Ifng", "Ccr5", "Ccr6", "Il21r", "Batf", "Irf4"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Il21", "Ifng", "Ccr5", "Ccr6", "Il21r", "Batf", "Irf4"), pt.size = 0.3)

# Tmem
VlnPlot(T_cell_subset, features = c("Ccr7", "Il21", "Icos", "Cd44"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Ccr7", "Il21", "Icos", "Cd44"), pt.size = 0.3)

# Treg
VlnPlot(T_cell_subset, features = c("Foxp3", "Ahr", "Lrrc32", "Il10", "Cd5"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Foxp3", "Ahr", "Lrrc32", "Il10", "Cd5"), pt.size = 0.3)


############ Trial and error #############

# Tfh
VlnPlot(T_cell_subset, features = c("Cxcr5", "Icos", "Cd84", "Il10", "Il21", "Il4", "Btla", "Cmaf"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Cxcr5", "Icos", "Cd84", "Il10", "Il21", "Il4", "Btla", "Cmaf"), pt.size = 0.3)

# Th22
VlnPlot(T_cell_subset, features = c("Ccr4", "Ccr6", "Ccr10", "Batf", "Stat3", "Il10", "Il21"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Ccr4", "Ccr6", "Ccr10", "Batf", "Stat3", "Il10", "Il21"), pt.size = 0.3)

# Th1
VlnPlot(T_cell_subset, features = c("Ccr5", "Cxcr3", "Ifngr2", "Il21r", "Il2", "Stat1", "Stat4"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Ccr5", "Cxcr3", "Ifngr2", "Il21r", "Il2", "Stat1", "Stat4"), pt.size = 0.3)

# Th2
VlnPlot(T_cell_subset, features = c("Ccr4", "Cxcr4", "Gata3", "Irf4", "Stat6", "Il4", "Il10", "Il21"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Ccr4", "Cxcr4", "Gata3", "Irf4", "Stat6", "Il4", "Il10", "Il21"), pt.size = 0.3)

# Th19
VlnPlot(T_cell_subset, features = c("Irf4", "Ccl17", "Ccl22", "Il10"), pt.size = 0.3, group.by = "orig.ident")
FeaturePlot(T_cell_subset, features = c("Irf4", "Ccl17", "Ccl22", "Il10"), pt.size = 0.3)

# DIFFERENTIAL GENE EXPRESSION -----------------------------------------------------------------
cd4_gex_073121 <- FindAllMarkers(T_cell_subset)
head(cd4_gex_073121)
write.csv(cd4_gex_073121, file= "/home/bmrc/Desktop/cd4_gex_073121.csv")

EnhancedVolcano(cd4_gex_073121,
                lab = rownames(cd4_gex_073121),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 10e-4,
                FCcutoff = 1.5,
                xlim = c(-5.5, 5.5),
                pointSize = 2.0,
                labSize = 3.0,
                title = 'Naive vs 14dpi',
                subtitle = NULL,
                caption = 'FC cutoff, 1.5; p-value cutoff, 10e-4',
                legendLabSize = 14,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
                colAlpha = 0.5,)