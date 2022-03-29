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

# 3. Annotation using SingleR package -----------------------------------------------------------------------------

integrated_annotation <- readRDS(file = "/home/phil/integrated_brucei/clustering_integrated.rds")
integrated_annotation

# load ImmGenData library
ref <- ImmGenData()
ref

"
# annotate integrated_annotation variable
whole_annotation <- SingleR(test = SummarizedExperiment(integrated_annotation), ref = ref, assay.type.test=1, 
			                                labels = ref$label.main)
whole_annotation
head(integrated_annotation@meta.data)

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
integrated_annotation$SingleR.pruned.calls <- whole_annotation$pruned.labels
integrated_annotation$SingleR.calls <- whole_annotation$labels

#Run UMAP------------------------------------------------------------------------------------------
integrated_annotation <- RunUMAP(integrated_annotation, dims = 1:10)
DimPlot(integrated_annotation, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.5)
head(integrated_annotation@meta.data)
integrated_annotation
"

#Annotation Diagnostics----------------------------------------------------------------------------
pred.grun <- SingleR(test=as.SingleCellExperiment(integrated_annotation), ref=ref, labels=ref$label.main, de.method="wilcox")
table(pred.grun$labels)
plotScoreHeatmap(pred.grun)
plotDeltaDistribution(pred.grun, ncol = 3)

all.markers <- metadata(pred.grun)$de.genes
sceG$labels <- pred.grun$labels


# Next, switch the identity class of all cells to reflect replicate ID
Idents(integrated_annotation) <- "SingleR.calls"

# How many cells are in each cluster
table(Idents(integrated_annotation))

# How many cells are in each cluster
integrated_annotation_counts <- table(Idents(integrated_annotation))
head(integrated_annotation_counts)

integrated_annotation_long <- as.data.frame(integrated_annotation_counts)
head(integrated_annotation_long)
colnames(integrated_annotation_long) <- c("cell_type", "Cell_count")
head(integrated_annotation_long)
rownames(integrated_annotation_long) <- integrated_annotation_long$x

integrated_annotation_long

bar_counts <- ggplot(integrated_annotation_long,            # Create ggplot2 plot scaled to 1.00
		                   aes(x = cell_type, y = Cell_count,
				                         fill = cell_type)) +
  geom_bar(stat = "identity", color ="black") + 
    #scale_y_continuous(labels = scales::percent_format()) + 
    scale_x_discrete(name ="Cell type") +
      theme(legend.position="right") +
        theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

bar_counts                               # Draw ggplot2 plot scaled to 1.00

# save annotated file
saveRDS(integrated_annotation, file = "annotation_integrated.rds")
