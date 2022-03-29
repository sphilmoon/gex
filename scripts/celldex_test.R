#Load the data
pyk <- readRDS(file = "/home/bmrc/Public/phil_ubuntu/sc/gex/pyk_brucei_4dpi/results_rds/harmony_integrated_111921_PM.rds")
pyk.diet <- DietSeurat(pyk)
pred.pyk <- as.SingleCellExperiment(pyk.diet)

# SingleR annotation using ImmGen reference dataset from the celldex package.
immgen <- ImmGenData()

whole_annotation <- SingleR(test = pred.pyk, 
                            ref = immgen, assay.type.test=1, 
                            labels = immgen$label.main)

pyk[["SingleR.labels"]] <- whole_annotation$labels

# Plot the annotated tables
DimPlot(pyk, reduction = "umap", group.by = "SingleR.labels", label = T, pt.size = 0.8)
FeaturePlot(pyk, features = c("Cd5l", "Cd36"))

##################### START #################################





########################## END ############################################
#Bcell annotation
Bcell_common_markers = c("Cd19", "Cd79b", "Ptprc", "Pax5")
FeaturePlot(combine, features = Bcell_common_markers) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
display.brewer.all() 