library(splatter)
library(Seurat)
params <- newSplatParams(batchCells = 2000, nGenes = 5000, 
                         group.prob = c(0.04,0.06,0.03,0.07,0.08,0.1,0.06,0.04,0.04,0.06,0.07,0.07,0.08,0.1,0.1), 
                         dropout.mid = 0.8,
                         dropout.shape = -0.8,
                         de.prob = 0.3,
                         de.facScale = c(0.11,0.12,0.11,0.11,0.11,0.15,0.12,0.11,0.13,0.12,0.15,0.13,0.13,0.12,0.11),
                         seed = 8)

sim <- splatSimulateGroups(params, verbose = FALSE)
simcounts <- counts(sim)
idents <- sim$Group

sim4 <- CreateSeuratObject(counts = counts(sim))
Idents(sim4) <- sim$Group

sim4 <- NormalizeData(sim4)
sim4 <- FindVariableFeatures(sim4, selection.method = "vst")
sim4 <- ScaleData(sim4)

sim4 <- RunPCA(sim4, features = VariableFeatures(object = sim4))
DimPlot(sim4, reduction = "pca")
sim4 <- RunTSNE(sim4, dims = 1:10)
DimPlot(sim4, reduction = "tsne")
sim4 <- RunUMAP(sim4, dims = 1:10)
DimPlot(sim4, reduction = "umap")

write.table(idents, file='idents3.csv', quote=FALSE, sep='\t', col.names = NA)

sim4 <- SCTransform(sim4, verbose = FALSE)
write.table(t(sim4@assays$SCT@scale.data), file='data3.csv', quote=FALSE, sep='\t', col.names = colnames(t(sim4@assays$SCT@scale.data)))

write.table(t(simcounts), file='data3_raw.csv', quote=FALSE, sep='\t', col.names = NA)
