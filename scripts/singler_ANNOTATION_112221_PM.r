#Load the data
pyk <- readRDS(file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/harmony_integrated_111921_PM.rds")
pyk

# SingleR annotation using ImmGen reference dataset from the celldex package.
immgen <- ImmGenData()
immgen

whole_annotation <- SingleR(test = SummarizedExperiment(pyk), 
                            ref = immgen, assay.type.test=1, 
                            labels = immgen$label.main)
whole_annotation
head(pyk@meta.data)

# Copy over the labels and pruned.labels 
# (Note: any other column of the results could be used as well)
pyk$SingleR.pruned.calls <- whole_annotation$pruned.labels
pyk$SingleR.calls <- whole_annotation$labels

# Save SingleR annotated file
#saveRDS(pyk, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/annotated_112221_PM.rds")

# Run UMAP
pyk <- RunUMAP(pyk, dims = 1:30)
DimPlot(pyk, reduction = "umap", group.by = "SingleR.calls", label = TRUE, pt.size = 0.5)

head(pyk@meta.data)
pyk

# Annotation Diagnostics
pred.grun <- SingleR(test=as.SingleCellExperiment(pyk), ref=immgen, 
					labels=immgen$label.main, de.method="wilcox")
table(pred.grun$labels)
plotScoreHeatmap(pred.grun)
plotDeltaDistribution(pred.grun, ncol = 3)

all.markers <- metadata(pred.grun)$de.genes
sceG$labels <- pred.grun$labels

# Next, switch the identity class of all cells to reflect replicate ID
Idents(pyk) <- "SingleR.calls"

# How many cells are in each cluster
table(Idents(pyk))

# How many cells are in each cluster
pyk_counts <- table(Idents(pyk))
head(pyk_counts)

pyk_long <- as.data.frame(pyk_counts)
head(pyk_long)
colnames(pyk_long) <- c("cell_type", "Cell_count")
head(pyk_long)
rownames(pyk_long) <- pyk_long$x
pyk_long

bar_pyk_long <- ggplot(pyk_long, # Create ggplot2 plot scaled to 1.00
		            aes(x = cell_type, y = Cell_count,
				    fill = cell_type)) + 
					geom_bar(stat = "identity", color ="black") + 
					# scale_y_continuous(labels = scales::percent_format()) + 
					scale_x_discrete(name ="Cell type") + 
					theme(legend.position="right") + 
					theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
bar_pyk_long # Draw ggplot2 plot scaled to 1.00
#### STOP HERE ####


# Can I create a Seurat object of just the Plasma cells and B cells?
B_Plasma_cells <- subset((pyk), idents = c("B cells", "B cells, pro"))

# How many cells are in each cluster
table(Idents(B_Plasma_cells))
B_Plasma_cells
head(B_Plasma_cells@meta.data)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(B_Plasma_cells) <- "seurat_clusters"
results <- SingleR(test = as.SingleCellExperiment(B_Plasma_cells), 
					ref = immgen, assay.type.test=1, 
					labels = immgen$label.fine) # OR I could use label.main
results

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
B_Plasma_cells$SingleR.pruned.calls_B <- results$pruned.labels
B_Plasma_cells$SingleR.calls_B <- results$labels
head(B_Plasma_cells@meta.data)

#Run UMAP
B_Plasma_cells <- RunUMAP(B_Plasma_cells, dims = 1:30)
DimPlot(B_Plasma_cells, reduction = "umap", label = TRUE)
DimPlot(B_Plasma_cells, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_Plasma_cells, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_Plasma_cells, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(B_Plasma_cells) <- "SingleR.calls_B"

# How many cells are in each cluster
table(Idents(B_Plasma_cells))

# Can I create a Seurat object of just the Plasma cells and B cells?
B_cells <- subset((B_Plasma_cells), 
				idents = c("B cells (B.CD19CONTROL)", "B cells (B.Fo)","B cells (B.FO)","B cells (B.FrE)", 
						"B cells (B.FRE)","B cells (B.FrF)","B cells (B.GC)","B cells (B.MEM)", 
						"B cells (B.MZ)","B cells (B.T1)","B cells (B.T2)","B cells (B.T3)", 
						"B cells (B1a)","B cells (B1A)","B cells (B1b)","B cells (preB.FrC)", 
						"B cells (preB.FrD)","B cells (preB.FRD)","B cells (proB.CLP)","B cells (proB.FrA)", 
						"B cells (proB.FrBC)"))
# How many cells are in each cluster
B_cell_counts <- table(Idents(B_cells), B_cells$group)
head(B_cell_counts)

B_cell_counts_long <- as.data.frame(B_cell_counts)    
colnames(B_cell_counts_long) <- c("B_cell", "Group", "Relative_abundance")
rownames(B_cell_counts_long) <- B_cell_counts_long$x
B_cell_counts_long

ggp <- ggplot(B_cell_counts_long, # Create ggplot2 plot scaled to 1.00 
				aes(x = Group, 
				y = Relative_abundance, 
				fill = B_cell)) + 
				geom_bar(position = "fill", stat = "identity", color ="black") + 
				scale_y_continuous(labels = scales::percent_format()) + 
				scale_x_discrete(name ="Group") + theme(legend.position="right")
ggp  # Draw ggplot2 plot scaled to 1.00

# Save counts in CSV
write.csv(B_cell_counts_long, file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/Bcell_counts_112221_PM.txt")

# Run UMAP
B_cells <- RunUMAP(B_cells, dims = 1:30)
DimPlot(B_cells, reduction = "umap", label = TRUE)
DimPlot(B_cells, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_cells, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(B_cells, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

# DOT PLOT MARKER GENES
cd_genes <- c("Cd19", "Cd93", "Ly6d", "Ebf1","Ms4a1","Vpreb3","Itgam","Zbtb32","Zbtb20", 
				"Cr2","Ighd","Fcer2a","Aicda","Mki67","Ighm","Sdc1","Ezh2","Bcl6", "Cd1d1", 
				"Cd9","Pax5","Ptprc", "Xbp1", "Jchain")
DotPlot(object = B_cells, features = cd_genes) + RotatedAxis() + 
		scale_colour_gradient2(low = "blue", mid = "white", high = "red")

# B Cells Only
FeaturePlot(B_cells, features = c("Cd79a", "Cd19", "Pax5", "Ebf1"), pt.size = 0.3) & 
			scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Cd79a", "Cd19", "Pax5", "Ebf1"), pt.size = 0.3, 
								group.by = "orig.ident")

# Transitional B cells
FeaturePlot(B_cells, features = c("Ly6d", "Ebf1", "Ms4a1", "Vpreb3"), pt.size = 0.3) & 
			scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Ly6d", "Ebf1", "Ms4a1", "Vpreb3"), pt.size = 0.3, 
								group.by = "orig.ident") 

# Transitional T2 B cells
FeaturePlot(B_cells, features = c("Ly6d", "Ms4a1", "Ighd", "Fcer2a"), pt.size = 0.3) & 
			scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Ly6d", "Ms4a1", "Ighd", "Fcer2a"), pt.size = 0.3, 
								group.by = "orig.ident") 

# Transitional T3 B cells
FeaturePlot(B_cells, features = c("Cd19","Cr2","Ighd","Fcer2a"), pt.size = 0.3) & 
			scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Cd19", "Cr2", "Ighd", "Fcer2a"), pt.size = 0.3, 
								group.by = "orig.ident")

# B1A B cells
FeaturePlot(B_cells, features = c("Zbtb20"), pt.size = 0.3) & 
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Zbtb20"), pt.size = 0.3, group.by = "orig.ident", 
								group.by = "orig.ident")

# B1B B cells
FeaturePlot(B_cells, features = c("Zbtb32", "Itgam", "Cd19"), pt.size = 0.3) & 
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Zbtb32", "Itgam", "Cd19"), pt.size = 0.3, 
								group.by = "orig.ident")

# Atypical Memory B cells
FeaturePlot(B_cells, features = c("Itgam","Zbtb32","Zbtb20", "Itgax", "Tbx21", "S100a6"), 
								pt.size = 0.3) &
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Itgam","Zbtb32","Zbtb20", "Itgax", "Tbx21", "S100a6"), 
								pt.size = 0.3, group.by = "orig.ident")

# Marginal Zone B cells
FeaturePlot(B_cells, features = c("Ly6d", "Ebf1", "Cr2", "Ms4a1", "Cd19"), pt.size = 0.3) &
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Ly6d", "Ebf1", "Cr2", "Ms4a1", "Cd19"), pt.size = 0.3, 
								group.by = "orig.ident")

# Follicular B cells
FeaturePlot(B_cells, features = c("Ighd", "Fcer2a", "Cd19", "Cr2"), pt.size = 0.3) & 
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Ighd", "Fcer2a", "Cd19", "Cr2"), pt.size = 0.3, 
								group.by = "orig.ident")

# Germinal center B cells
FeaturePlot(B_cells, features = c("Aicda", "Mki67"), pt.size = 0.3) & 
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Aicda", "Mki67"), pt.size = 0.3, 
								group.by = "orig.ident")

# Plasma B cells
FeaturePlot(B_cells, features = c("Sdc1", "Xbp1", "Ighm", "Jchain"), pt.size = 0.3) & 
								scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(B_cells, features = c("Sdc1", "Xbp1", "Ighm", "Jchain"), pt.size = 0.3, 
								group.by = "orig.ident")