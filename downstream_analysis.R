library(ggplot2)
library(dplyr)
library(Seurat)
library(patchwork)

latent <- read.table("IPMN_2d.csv", sep=",", header=TRUE)
p <- ggplot(data = latent) + geom_point(mapping = aes(x = X0, y = X1))
p

#latent <- read.table("IPMN_encoded.csv", sep=",", header=TRUE)
expression <- read.table("GSM5032771_ZX-2.expression_matrix.txt", sep="\t", header=TRUE)

selected.data <- seq(from = 1, to = nrow(latent), by = 2)
latent <- latent[selected.data,]
expression <- expression[,selected.data]
rownames(latent) <- colnames(expression)
colnames(latent) <- paste0("SAE_", 1:2)

ipmn <- CreateSeuratObject(counts = expression, project = "IPMN", min.cells = 3)
ipmn
VlnPlot(ipmn, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(ipmn, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

ipmn <- NormalizeData(ipmn, normalization.method = "LogNormalize", scale.factor = 10000)

ipmn <- FindVariableFeatures(ipmn, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(ipmn), 10)
plot1 <- VariableFeaturePlot(ipmn)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(ipmn)
ipmn <- ScaleData(ipmn, features = all.genes)
#ipmn <- SCTransform(ipmn, verbose = FALSE)

ipmn <- RunPCA(ipmn, features = VariableFeatures(object = ipmn))

ipmn <- FindNeighbors(ipmn, dims = 1:10)
ipmn <- FindClusters(ipmn, resolution = 0.5)

sae <- CreateDimReducObject(embeddings = as.matrix(latent), key = "SAE_", assay = DefaultAssay(ipmn))
ipmn[["sae"]] <- sae
DimPlot(ipmn, reduction = "sae", pt.size = 0.5)

ipmn <- ProjectDim(ipmn, reduction = "sae")
DimHeatmap(ipmn, reduction = "sae", dims = 1, cells = 500, projected = TRUE, balanced = TRUE)

ipmn.markers <- FindAllMarkers(ipmn, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
select.markers <- ipmn.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

VlnPlot(ipmn, features = c("IL1B", "LYZ"))
VlnPlot(ipmn, features = c("IGKC", "IGHG2"))
VlnPlot(ipmn, features = c("MS4A1", "CD37"))
VlnPlot(ipmn, features = c("CD3D", "CD2"))
VlnPlot(ipmn, features = c("COL1A1", "LUM"))
VlnPlot(ipmn, features = c("PI3"))
VlnPlot(ipmn, features = c("SPARCL1", "CLDN5"))
VlnPlot(ipmn, features = c("RGS5", "REG3A"))
VlnPlot(ipmn, features = c("LIF"))

feature.gene <- c("IL1B", "LYZ", "IGKC", "IGHG2", "MS4A1", "CD37", "CD3D", "CD2", "COL1A1", "LUM", "PI3", "SPARCL1", "CLDN5", "RGS5", "REG3A", "CD14", "FCGR3A", "MS4A7")
FeaturePlot(ipmn, features = feature.gene)
DotPlot(ipmn, features = feature.gene)

new.cluster.ids <- c("Monocytes", "Fibroblasts", "Plasma cells", "Monocytes", "T cells", "Macrophage", "B cells", "B cells", "Endothelial cells", "Plasma cells")
names(new.cluster.ids) <- levels(ipmn)
ipmn <- RenameIdents(ipmn, new.cluster.ids)
DimPlot(ipmn, reduction = "sae", label = TRUE, pt.size = 0.5) + NoLegend()


