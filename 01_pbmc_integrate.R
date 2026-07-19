setRepositories(ind = 1:3)
required.packages <- c("Seurat", "BiocManager", "harmony", "glmGamPoi", "devtools")
InstallMissing <- function(packages) {
    for (pkg in packages) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, dependencies = TRUE)
            library(pkg, character.only = TRUE)
        }
    }
}

BiocManager::install("batchelor")
library(Seurat)
library(dplyr)
devtools::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)


pbmc3k <- readRDS("data/pbmc_integration/pbmc3k_annotated.rds")
pbmc3k$dataset <- "pbmc3k"
pbmc4k <- readRDS("data/pbmc_integration/pbmc4k.rds")
pbmc4k$dataset <- "pbmc4k"

# Merge the two datasets
merged <- merge(pbmc3k, y = pbmc4k, add.cell.ids = c("pbmc3k", "pbmc4k"))

# Perform standard workflow for visualization and clustering
merged <- NormalizeData(merged) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(npcs = 30) %>%
  FindNeighbors(dims = 1:20) %>%
  RunUMAP(dims = 1:20)

# Visualize the merged dataset
# do the datasets merge into each other well?
DimPlot(merged, reduction = "umap", group.by = "dataset")

merged <- FindClusters(merged, resolution = 0.3)
merged$merged_cluster <- Idents(merged)

# visualise the clusters
DimPlot(merged, reduction = "umap", group.by = "merged_cluster", label = TRUE)

DimPlot(merged, reduction = "umap", group.by = c("dataset",
                                                 "merged_cluster"), label = TRUE)


# Let's integrate
merged <- JoinLayers(merged)
merged[["RNA"]] <- split(merged[["RNA"]], f = merged$dataset)
merged <- NormalizeData(ifnb) %>% FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

integrated <-  IntegrateLayers(object = merged, method = FastMNNIntegration,
                                   orig.reduction = "pca",
                                   new.reduction = "integrated.mnn",
                                          verbose = FALSE)
integrated <- IntegrateLayers(
    object = integrated, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "integrated.harmony",
    verbose = FALSE
)


integrated <- IntegrateLayers(
    object = integrated, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE
)

integrated <- FindNeighbors(integrated, reduction = "integrated.cca", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.3, cluster.name = "cca_clusters")
integrated <- RunUMAP(integrated, reduction = "integrated.cca", dims = 1:20,
                      reduction.name = "umap.cca")

DimPlot(integrated, reduction = "umap.cca", group.by = c("dataset",
                                                 "merged_cluster", "cca_clusters"),
        label = TRUE)

integrated <- FindNeighbors(integrated, reduction = "integrated.harmony", dims = 1:20)
integrated <- FindClusters(integrated, resolution = 0.3,
                           cluster.name = "harmony_clusters")
integrated <- RunUMAP(integrated, reduction = "integrated.harmony", dims = 1:20,
                      reduction.name = "umap.harmony")

DimPlot(integrated, reduction = "umap.harmony", group.by = c("dataset",
                                                         "harmony_clusters", "cca_clusters"),
        label = TRUE)

# We can use the same idea to 'transfer' annotation from pbmc3k to pbmc4k

pbmc3k <- NormalizeData(pbmc3k) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(dims = 1:20) %>%
    RunUMAP(dims = 1:20)
pbmc4k <- NormalizeData(pbmc4k) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData()  %>%
    RunPCA(npcs = 30) %>%
    FindNeighbors(dims = 1:20) %>%
    RunUMAP(dims = 1:20)

anchors <- FindTransferAnchors(reference = pbmc3k,
                               query = pbmc4k,
                               dims = 1:30,
                               reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = pbmc3k$seurat_annotations,
                            dims = 1:30)
pbmc4k <- AddMetaData(pbmc4k, metadata = predictions)
DimPlot(pbmc4k, group.by = "predicted.id")

# uamp is only for visualisatio
# always check using markers if the results make sense


markers.B <- FindMarkers(object = pbmc4k, ident.1 = "B", only.pos=T, group.by = "predicted.id")

