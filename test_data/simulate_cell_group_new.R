library(splatter)
library(Seurat)
params <- newSplatParams(batchCells = 2000, nGenes = 5000, 
                         group.prob = c(0.09,0.09,0.1,0.08,0.12,0.11,0.1,0.1,0.11,0.1),
                         #dropout.type = "experiment",
                         dropout.mid = 0.8,
                         dropout.shape = -0.8,
                         de.prob = 0.3,
                         de.facScale = c(0.1,0.12,0.11,0.13,0.11,0.09,0.07,0.1,0.11,0.08),
                         seed = 8)

sim <- splatSimulateGroups(params, verbose = FALSE)
simcounts <- counts(sim)
idents <- sim$Group

sim4 <- CreateSeuratObject(counts = counts(sim))
Idents(sim4) <- sim$Group

sim4 <- NormalizeData(sim4)
sim4 <- FindVariableFeatures(sim4, selection.method = "vst")
sim4 <- ScaleData(sim4)

sim4 <- RunPCA(sim4)
DimPlot(sim4, reduction = "pca")
sim4 <- RunTSNE(sim4, dims = 1:10)
DimPlot(sim4, reduction = "tsne")
sim4 <- RunUMAP(sim4, dims = 1:20)
DimPlot(sim4, reduction = "umap")

write.table(idents, file='data3_dropout_idents_0.csv', quote=FALSE, sep='\t', col.names = NA)

sim4 <- SCTransform(sim4, verbose = FALSE)
write.table(sim4@assays$SCT@scale.data, file='data_SCT.csv', quote=FALSE, sep='\t', col.names = TRUE, row.names = TRUE)

write.table(t(simcounts), file='data_raw.csv', quote=FALSE, sep='\t', col.names = NA)

write.table(simcounts, file='data_raw_t.csv', quote=FALSE, sep=',', col.names = TRUE, row.names = TRUE)






