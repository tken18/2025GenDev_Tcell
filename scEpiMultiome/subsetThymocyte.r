#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO_DATA/integrated.data.rds')
subsetData <- subset(integrated.data, idents = c('0', '1'))


# RNA -----
rna <- CreateSeuratObject(subsetData@assays$RNA)
rna <- AddMetaData(rna, subsetData@meta.data)

sample.list <- SplitObject(rna, split.by = 'genotype')
for (i in 1:length(sample.list)) {
  sample.list[[i]] <- SCTransform(sample.list[[i]], verbose = FALSE)
  
}

features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 2000)
list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = list, normalization.method = 'SCT', anchor.features = features)
integrated.rna <- IntegrateData(anchorset = anchors, normalization.method = 'SCT')
integrated.rna <- RunPCA(integrated.rna, verbose = FALSE)
integrated.rna <- RunUMAP(integrated.rna, reduction = 'pca', dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# ATAC -----
atac <- CreateSeuratObject(subsetData@assays$ATAC)
atac <- AddMetaData(atac, subsetData@meta.data)

sample.list <- SplitObject(atac, split.by = 'genotype')
for (i in 1:length(sample.list)) {

  sample.list[[i]] <- FindTopFeatures(sample.list[[i]], min.cutoff = 'q0')
  sample.list[[i]] <- RunTFIDF(sample.list[[i]])
  sample.list[[i]] <- RunSVD(sample.list[[i]])
  
}

combined.atac <- merge(sample.list[[1]], sample.list[[2]])
combined.atac <- FindTopFeatures(combined.atac, min.cutoff = 10)
combined.atac <- RunTFIDF(combined.atac)
combined.atac <- RunSVD(combined.atac)
combined.atac <- RunUMAP(combined.atac, reduction = 'lsi', dims = 2:30)

integration.anchors <- FindIntegrationAnchors(object.list = sample.list, anchor.features = 2000, reduction = 'rlsi', dims = 2:30)
integrated.atac <- IntegrateEmbeddings(anchorset = integration.anchors, reductions = combined.atac[['lsi']], new.reduction.name = 'integrated_lsi', k.weight = 10, dims.to.integrate = 1:30)
integrated.atac <- RunUMAP(integrated.atac, reduction = 'integrated_lsi', dims = 2:30, reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')


# WNN -----
data <- integrated.rna
data[['ATAC']] <- integrated.atac[['RNA']]
data@reductions$lsi <- integrated.atac@reductions$integrated_lsi
data@reductions$umap.atac <- integrated.atac@reductions$umap.atac

data <- FindMultiModalNeighbors(data, reduction.list = list('pca', 'lsi'), dims.list = list(1:50, 2:50))
data <- RunUMAP(data, nn.name = 'weighted.nn', reduction.name = 'wnn.umap', reduction.key = 'wnnUMAP_')
data <- FindClusters(data, graph.name = 'wsnn', algorithm = 3, verbose = FALSE, resolution = 0.5)
data$genotype <- factor(data$genotype, levels = c('WT', 'dKO'))


# Run chromVAR -----
DefaultAssay(data) <- 'ATAC'
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
keep.peaks <- as.logical(seqnames(granges(data)) %in% main.chroms)
data[["ATAC"]] <- subset(data[["ATAC"]], features = rownames(data[["ATAC"]])[keep.peaks])
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
data <- AddMotifs(object = data, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

saveRDS(data, 'thymocyteWNN.rds')


# Peak calling -----
DefaultAssay(data) <- 'ATAC'
peaks <- CallPeaks(object = data, group.by = 'genotype', idents = c(), macs2.path = '/usr/local/bin/macs2')
write.table(peaks, 'WholeCellPeaks.txt', row.names = T, col.names = T, quote = F, sep = '\t')

Idents(data) <- 'genotype'
peaks <- CallPeaks(object = data, group.by = 'seurat_clusters', idents = 'WT', macs2.path = '/usr/local/bin/macs2')
write.table(peaks, 'WTpeaks.txt', row.names = T, col.names = T, quote = F, sep = '\t')

peaks <- CallPeaks(object = data, group.by = 'seurat_clusters', idents = 'dKO', macs2.path = '/usr/local/bin/macs2')
write.table(peaks, 'dKOpeaks.txt', row.names = T, col.names = T, quote = F, sep = '\t')


# Stacking bar plots -----
pal.stack$cluster <- as.character(0:11)
pal.stack$cluster <- factor(pal.stack$cluster, levels = c('11', '4', '0', '1', '2', '7', '5', '6', '9', '3', '10', '8', '12'), ordered = TRUE)
pal.stack <- pal.stack[order(pal.stack$cluster),]

ratio <- as.data.frame.matrix(prop.table(table(Idents(data), data$genotype), margin = 2))
ratio_g = gather(ratio, key = group, value = ratio)
ratio_g$cluster <- levels(data@active.ident)
ratio_g$cluster <- factor(ratio_g$cluster, levels = levels(data@active.ident))
ratio_g$cluster <- factor(ratio_g$cluster, levels = c('11', '4', '0', '1', '2', '7', '5', '6', '9', '3', '10', '8'))
ratio_g$group = factor(ratio_g$group, levels = c('WT', 'dKO'))

p = NULL
p <- ggplot(data = ratio_g, aes(x = group, y = ratio, fill = cluster, label = cluster))
p <- p + geom_bar(stat = 'identity', position = 'fill', col = 'black', width = 0.7)
p <- p + theme_minimal() + coord_fixed(ratio = 3)
p

