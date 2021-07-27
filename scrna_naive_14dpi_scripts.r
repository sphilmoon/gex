library(dplyr)
library(Seurat)
library(patchwork)
library(SoupX)
library(ggplot2)

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


### INSERT JUPYTER NOTEBOOK HERE ### 




# 3. Annotation using SingleR package -----------------------------------------------------------------------------


infected_umap <- readRDS(file = "/home/phil/infected/seurat/14dpi_umap_trial.rds")
infected_umap

# load ImmGenData library
ref <- ImmGenData()
ref

# annotate infected_umap file
infected_annotation <- SingleR(test = SummarizedExperiment(infected_umap), ref = ref, assay.type.test=1, 
			                                   labels = ref$label.main)
infected_annotation
head(infected_umap@meta.data)

# Copy over the labels and pruned.labels (Note: any other column of the results could be used as well)
infected_umap$SingleR.pruned.calls <- infected_annotation$pruned.labels
infected_umap$SingleR.calls <- infected_annotation$labels











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

infected_data <- CreateSeuratObject(counts = infected_soupx, project = "infected", min.cells = 3, min.features = 200)

