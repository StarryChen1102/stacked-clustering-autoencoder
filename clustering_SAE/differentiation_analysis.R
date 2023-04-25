library(Seurat)
library(monocle3)
library(ggplot2)

latent <- read.table("bm4_diff/bm4subset_model_pre_umap.csv", sep=",", header=T)
p <- ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1))
p

id <- read.table("bm4_diff/bm4subset_idents.csv", sep="\t", header=TRUE)
expression <- read.table("bm4_diff/bm4subset_raw.csv", sep=",", header=TRUE)
#rownames(expression) <- expression$X
#expression <- expression[2:19568]
#expression <- as.data.frame(t(expression))

rownames(latent) <- colnames(expression)
colnames(latent) <- paste0("SAE_", 1:2)

bm4 <- CreateSeuratObject(counts = expression, project = "data1", min.cells = 3)
bm4
VlnPlot(bm4, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(bm4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Idents(bm4) <- id$x

bm4 <- SCTransform(bm4)
bm4 <- RunPCA(bm4, features = VariableFeatures(object = bm4))
DimPlot(bm4, reduction = "pca", pt.size = 0.5)
bm4 <- RunUMAP(bm4, dims = 1:30)
DimPlot(bm4, reduction = "umap")
bm4 <- RunTSNE(bm4, dims = 1:30)
DimPlot(bm4, reduction = "tsne")

sae <- CreateDimReducObject(embeddings = as.matrix(latent), key = "SAE_", assay = DefaultAssay(bm4))
bm4[["sae"]] <- sae
DimPlot(bm4, reduction = "sae", pt.size = 0.5)

bm4 <- ProjectDim(bm4, reduction = "sae")
DimHeatmap(bm4, reduction = "sae", dims = 1, cells = 500, projected = TRUE, balanced = TRUE)

data <- GetAssayData(bm4, assay = 'RNA', slot = 'counts')
cell_metadata <- bm4@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method="tSNE")
cds <- reduce_dimension(cds, reduction_method="UMAP")

#cds.embed <- cds@int_colData$reducedDims$tSNE
int.embed <- Embeddings(bm4, reduction = "sae")
#int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds@int_colData$reducedDims$tSNE <- int.embed
cds <- cluster_cells(cds, reduction_method = "UMAP") 

cds <- learn_graph(cds)
head(colData(cds))

cds = order_cells(cds) 

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           group_label_size=4,
           cell_size=1)  

plot_cells(cds)

plot_genes_in_pseudotime(cds["S100A9"], color_cells_by = "pseudotime", min_expr = 0.5)







