#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO_DIR/thymocyteWNN.rds')
WT <- subset(data, genotype == 'WT')
dKO <- subset(data, genotype == 'dKO')

threshold <- 0.5 # for Cicero


# Run Monocle3 using whole cells -----
cds <- as.cell_data_set(subset(data, idents = c('8', '6', '5'), invert = T), assay = 'ATAC')
reducedDim(cds, type = 'UMAP') <- cds@int_colData@listData$reducedDims$WNN.UMAP
cds <- cluster_cells(cds, reduction_method = 'UMAP', resolution = 0.001)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = which.max(data@assays$SCT@scale.data['Flt3',]))
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3.5)
data <- AddMetaData(data, cds@principal_graph_aux@listData[['UMAP']][['pseudotime']], 'monocle')
cds_subset <- choose_cells(cds) #Select transient cells using a graphical user interface according to the Monocle3 tutorials.

saveRDS(cds_subset, 'subsetTransientCells.rds')
saveRDS(data, 'thymocyteWNN_monocle3.rds')


# Run Monocle3 and Cicero using WT cells -----
cds <- as.cell_data_set(subset(WT, idents = c('8', '6', '5'), invert = T), assay = 'ATAC')
reducedDim(cds, type = 'UMAP') <- cds@int_colData@listData$reducedDims$WNN.UMAP
cds <- cluster_cells(cds, reduction_method = 'UMAP', resolution = 0.004)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = which.max(WT@assays$SCT@scale.data['Flt3',]))
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3.5)
WT <- AddMetaData(WT, cds@principal_graph_aux@listData[['UMAP']][['pseudotime']], 'monocle')
saveRDS(WT, 'thymocyteWNN_monocle3WT.rds')

cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$WNN.UMAP)
DefaultAssay(WT) <- 'ATAC'
genome <- seqlengths(WT)
genome.df <- data.frame('chr' = names(genome), 'length' = genome)
conns <- run_cicero(cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans, threshold = threshold)
Links(WT) <- links
saveRDS(WT, 'cicero05WT.rds')


# Run Monocle3 and Cicero using dKO cells -----
cds <- as.cell_data_set(subset(dKO, idents = c('8', '6', '5'), invert = T), assay = 'ATAC')
reducedDim(cds, type = 'UMAP') <- cds@int_colData@listData$reducedDims$WNN.UMAP
cds <- cluster_cells(cds, reduction_method = 'UMAP', resolution = 0.004)
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = which.max(dKO@assays$SCT@scale.data['Flt3',]))
plot_cells(cds, color_cells_by = 'pseudotime', label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 3.5)
dKO <- AddMetaData(dKO, cds@principal_graph_aux@listData[['UMAP']][['pseudotime']], 'monocle')
saveRDS(dKO, 'thymocyteWNN_monocle3dKO.rds')

cicero <- make_cicero_cds(cds, reduced_coordinates = reducedDims(cds)$WNN.UMAP)
DefaultAssay(dKO) <- 'ATAC'
genome <- seqlengths(dKO)
genome.df <- data.frame('chr' = names(genome), 'length' = genome)
conns <- run_cicero(cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans, threshold = threshold)
Links(dKO) <- links
saveRDS(dKO, 'cicero05dKO.rds')


## Rag ##
range <- c('chr2-101550998-101552999', 'chr2-101598685-101600974', 'chr2-101665094-101665827', 'chr2-101649348-101650864')
range.show = StringToGRanges(regions = range)
CoveragePlot(WT, region = 'chr2-101529844-101699843', assay = 'ATAC', peaks = F, ranges = range.show, ymax = 130) & scale_fill_manual(values = pal2)
CoveragePlot(dKO, region = 'chr2-101529844-101699843', assay = 'ATAC', peaks = F, ranges = range.show, ymax = 130)

## Id2 ##
CoveragePlot(WT, region = 'chr12-25070000-25226200', assay = 'ATAC', peaks = F, ymax = 220)
CoveragePlot(dKO, region = 'chr12-25070000-25226200', assay = 'ATAC', peaks = F, ymax = 220)


