library(ggplot2)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(Seurat)
library(scran)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

# Load annotated rds file
annotation_integrated <- readRDS(file = "/home/phil/integrated_brucei/annotation_integrated.rds")

# 1st, filter: Subset T cells -----------------------------------------------
# Can I create a Seurat object of just the Plasma cells and B cells?
T_cells_1st <- subset((annotation_integrated), idents = "T cells")

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
																		                                                        "T cells (T.CD4.24H)", "T cells (T.4Nve)", "T cells (T.4NVE)", 
																									                                                  "T cells (T.CD4.CTR)", "T cells (T.CD4CONTROL)",
																									                                                  "T cells (T.4.PLN)", "T cells (T.4.Pa)", "T cells (T.4FP3+25+)"))
                        
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
                        
#Save counts in CSV to see my subset annotation
#write.csv(T_cell_counts, file = "/home/Public/cd4_subset_annotation")

#Run UMAP
T_cells_umap <- RunUMAP(T_cell_subset, dims = 1:10)
DimPlot(T_cell_subset, reduction = "umap", label = TRUE)
DimPlot(T_cell_subset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cell_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(T_cell_subset, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
T_cell_markers <- FindAllMarkers(T_cell_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
T_cell_markers %>%
	  group_by(cluster) %>%
	    top_n(n = 10, wt = avg_log2FC)

    head(n=100, T_cell_markers)

    # T.4.PLN T cells genes
    VlnPlot(T_cell_subset, features = c("Ccr9", "Cd8b1", "Cd8a", "Itgae", "Rgcc", "Klrd1", "Fam241a"), pt.size = 0.3, group.by = "orig.ident")  
    FeaturePlot(T_cell_subset, features = c("Ccr9", "Cd8b1", "Cd8a", "Itgae", "Rgcc", "Klrd1", "Fam241a"), pt.size = 0.3) & 
	      scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

      # look at the most gex in my data, then compare to the literature
      # Next, plot the most expressed genes in the CD4 subset.

