#!/bin/Rscript

setwd('/WD_PATH/')

# Input data -----
WT.rna <- readRDS('/PATH_TO_DIR/WTfilteredRNA.rds')
WT.atac <- readRDS('/PATH_TO_DIR/WTfilteredATAC.rds')
dKO.rna <- readRDS('/PATH_TO_DIR/dKOfilteredRNA.rds')
dKO.atac <- readRDS('/PATH_TO_DIR/dKOfilteredATAC.rds')


# RNA data processing -----
WT.rna <- SCTransform(WT.rna, verbose = FALSE)
dKO.rna <- SCTransform(dKO.rna, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = c(WT.rna, dKO.rna), nfeatures = 3000)
list <- PrepSCTIntegration(object.list = c(WT.rna, dKO.rna), anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = 'SCT', anchor.features = features)
integrated.rna <- IntegrateData(anchorset = anchors, normalization.method = 'SCT')
integrated.rna <- RunPCA(integrated.rna, verbose = FALSE)
integrated.rna <- RunUMAP(integrated.rna, reduction = 'pca', dims = 1:30)


# ATAC data processing -----
# We will only use peaks in standard chromosomes
## Normalization ##
WT.atac <- FindTopFeatures(WT.atac, min.cutoff = 'q0')
WT.atac <- RunTFIDF(WT.atac)
WT.atac <- RunSVD(WT.atac)

dKO.atac <- FindTopFeatures(dKO.atac, min.cutoff = 'q0')
dKO.atac <- RunTFIDF(dKO.atac)
dKO.atac <- RunSVD(dKO.atac)

## Integration ##
combined.atac <- merge(WT.atac, dKO.atac)
combined.atac <- FindTopFeatures(combined.atac, min.cutoff = 10)
combined.atac <- RunTFIDF(combined.atac)
combined.atac <- RunSVD(combined.atac)
combined.atac <- RunUMAP(combined.atac, reduction = 'lsi', dims = 2:30)

hm.integrated_all <- RunHarmony(
  object = combined.atac,
  group.by.vars = 'genotype',
  reduction = 'lsi',
  asay.use = 'peaks',
  project.dim = FALSE
)

DefaultAssay(hm.integrated_all) <- 'RNA'
hm.integrated_all <- RunUMAP(object = hm.integrated_all, dims = 1:15, reduction = 'harmony')
hm.integrated_all <- FindNeighbors(object = hm.integrated_all,reduction = 'harmony', dims = 1:20)
hm.integrated_all <- FindClusters(object　=hm.integrated_all,resolution = 0.5)


# WNN -----
integrated.data[['ATAC']] <- integrated.atac[['RNA']]
integrated.data@reductions$lsi <- integrated.atac@reductions$lsi
integrated.data@reductions$harmony <- integrated.atac@reductions$harmony

## RNA ##
DefaultAssay(integrated.data) <- 'RNA'
integrated.data <- RunUMAP(integrated.data, dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

## ATAC ##
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(integrated.data) <- 'ATAC'
integrated.data <- RunUMAP(integrated.data, reduction = 'harmony', dims = 1:15, reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')

integrated.data <- FindMultiModalNeighbors(integrated.data, reduction.list = list('pca', 'harmony'), dims.list = list(1:50, 2:50))
integrated.data <- RunUMAP(integrated.data, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')
integrated.data <- FindClusters(integrated.data, graph.name = 'wsnn', algorithm = 3, verbose = FALSE, resolution = 0.2)
integrated.data$genotype <- factor(integrated.data$genotype, levels = c('WT', 'dKO'))
saveRDS(integrated.data, 'integrated.data.rds')


p1 <- DimPlot(integrated.data, reduction = 'umap.rna', label.size = 5, group.by = 'genotype', cols = c('#3B9AB2', '#F21A00')) + ggtitle('RNA') + NoLegend()
p2 <- DimPlot(integrated.data, reduction = 'umap.atac', label.size = 5, group.by = 'genotype', cols = c('#3B9AB2', '#F21A00')) + ggtitle('ATAC')
p3 <- DimPlot(integrated.data, reduction = 'wnn.umap', label.size = 5, group.by = 'genotype', cols = c('#3B9AB2', '#F21A00')) + ggtitle('WNN') + NoLegend()
p3 + p1 + p2 & theme(plot.title = element_text(hjust = 0.5))

p1 <- DimPlot(integrated.data, reduction = 'umap.rna', label.size = 5, repel = TRUE, label = TRUE, split.by = 'genotype') + ggtitle('RNA') + NoLegend()
p2 <- DimPlot(integrated.data, reduction = 'umap.atac', label.size = 5, repel = TRUE, label = TRUE, split.by = 'genotype') + ggtitle('ATAC')
p3 <- DimPlot(integrated.data, reduction = 'wnn.umap', label.size = 5, repel = TRUE, label = TRUE, split.by = 'genotype') + ggtitle('WNN') + NoLegend()
p3 + p1 + p2 & theme(plot.title = element_text(hjust = 0.5))

p1 <- DimPlot(integrated.data, reduction = 'umap.rna', label.size = 5, repel = TRUE, label = TRUE) + ggtitle('RNA') + NoLegend()
p2 <- DimPlot(integrated.data, reduction = 'umap.atac', label.size = 5, repel = TRUE, label = TRUE) + ggtitle('ATAC')
p3 <- DimPlot(integrated.data, reduction = 'wnn.umap', label.size = 5, repel = TRUE, label = TRUE) + ggtitle('WNN') + NoLegend()
p3 + p1 + p2 & theme(plot.title = element_text(hjust = 0.5))
