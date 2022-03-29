# Load the data
pyk_integrated <- readRDS(file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/harmony_integrated_111921_PM.rds")
pyk_integrated

annot <- ImmGenData()
annot

# Annotation Diagnostics
pred.grun <- SingleR(test=as.SingleCellExperiment(pyk_integrated), ref=annot, assay.type.test=1, 
                     labels=annot$label.main, de.method="wilcox")
table(pred.grun$labels)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(pyk_integrated) <- "SingleR.calls"

# How many cells are in each cluster
table(Idents(pyk_integrated))

# How many cells are in each cluster
pyk_integrated_counts <- table(Idents(pyk_integrated))
head(pyk_integrated_counts)

# to visualize the annotated table
pyk_integrated_long <- as.data.frame(pyk_integrated_counts)
head(pyk_integrated_long)
colnames(pyk_integrated_long) <- c("cell_type", "Cell_count")
head(pyk_integrated_long)
rownames(pyk_integrated_long) <- pyk_integrated_long$x
pyk_integrated_long

# Create ggplot2 plot scaled to 1.00
bar_counts <- ggplot(pyk_integrated_long, aes(x = cell_type, y = Cell_count, fill = cell_type)) + 
  geom_bar(stat = "identity", color ="black") + scale_x_discrete(name ="Cell type") + 
  theme(legend.position="right") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Draw ggplot2 plot scaled to 1.00
bar_counts

# save annotated file
#saveRDS(pyk_integrated, file = "pyk_annotation_integrated_112121_PM.rds")


######################### WORKS WELL UNTIL HERE #########################


# Can I create a Seurat object of just the Plasma cells and B cells?
b_cells_general <- subset((pyk_integrated), idents = c("B cells", "B cells, pro"))

# How many cells are in each cluster
table(Idents(b_cells_general))
b_cells_general
head(b_cells_general@meta.data)

# Next, switch the identity class of all cells to reflect replicate ID
Idents(b_cells_general) <- "seurat_clusters"
results <- SingleR(test = as.SingleCellExperiment(b_cells_general), 
                   ref = annot, assay.type.test=1, 
                   labels = annot$label.fine) # OR I could use label.main
results

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
b_cells_general$SingleR.pruned.calls_B <- results$pruned.labels
b_cells_general$SingleR.calls_B <- results$labels
head(b_cells_general@meta.data)

#Run UMAP
b_cells_general <- RunUMAP(b_cells_general, dims = 1:30)
DimPlot(b_cells_general, reduction = "umap", label = TRUE)
DimPlot(b_cells_general, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_general, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_general, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)



#Re-run UMAP on new created Bcell object
b_cells_general <- NormalizeData(b_cells_general) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
ElbowPlot(b_cells_general, n=50)
b_cells_general <- b_cells_general %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)
b_cells_general <- b_cells_general %>% 
  RunUMAP(reduction = "harmony", dims = 1:25) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:25) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()
DimPlot(b_cells_general, reduction = "umap",  label = T, pt.size = 1)


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
# How many cells are in each cluster
B_cell_counts <- table(Idents(b_cells_subset), b_cells_subset$group)
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
b_cells_subset <- RunUMAP(b_cells_subset, dims = 1:30)
DimPlot(b_cells_subset, reduction = "umap", label = TRUE)
DimPlot(b_cells_subset, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_subset, reduction = "umap", group.by = "orig.ident", label = TRUE, pt.size = 0.5)
DimPlot(b_cells_subset, reduction = "umap", group.by = "SingleR.calls_B", label = TRUE, pt.size = 0.5)

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
FeaturePlot(b_cells_subset, features = c("Zbtb20"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Zbtb20"), pt.size = 0.3, group.by = "orig.ident", 
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
FeaturePlot(b_cells_subset, features = c("Ighd", "Fcer2a", "Cd19", "Cr2"), pt.size = 0.3) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
VlnPlot(b_cells_subset, features = c("Ighd", "Fcer2a", "Cd19", "Cr2"), pt.size = 0.3, 
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
