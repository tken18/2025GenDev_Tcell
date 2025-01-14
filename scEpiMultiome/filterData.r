#!/bin/Rscript

setwd('/WD_PATH/')
dataDir <- '/PATH_TO_H5_DIRECTORY/'
genotype <- c('WT', 'dKO')

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- 'mm10'

for (i in 1:length(genotype)){

input <- Read10X_h5(paste0(dataDir, genotype[i], '/filtered_feature_bc_matrix.h5'))

rna_counts <- input$`Gene Expression`
atac_counts <- input$Peaks

data <- CreateSeuratObject(counts = rna_counts)
data[['percent.mt']] <- PercentageFeatureSet(data, pattern = '^mt-')
data$genotype <- genotype[i]

grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]

frag.file <- paste0(dataDir, genotype[i], '/atac_fragments.tsv.gz')
chrom_assay <- CreateChromatinAssay(counts = atac_counts, sep = c(':', '-'), genome = 'mm10', fragments = frag.file, min.cells = 10, annotation = annotations)

data[['ATAC']] <- chrom_assay
VlnPlot(rna, features = c('nCount_ATAC', 'nFeature_RNA', 'nCount_RNA','percent.mt'), ncol = 3, log = TRUE, pt.size = 0.05) + NoLegend()

filtered <- subset(x = data, subset = nCount_ATAC < 1e5 & nCount_ATAC > 5e3 & nCount_RNA < 35000 & nCount_RNA > 3500 & nFeature_RNA < 8000 & nFeature_RNA > 500 & percent.mt < 25)


filtered.rna <- filtered@assays$RNA
filtered.rna <- CreateSeuratObject(filtered.rna)
filtered.rna@meta.data <- filtered@meta.data

filtered.atac <- filtered@assays$ATAC
filtered.atac <- CreateSeuratObject(filtered.atac)
filtered.atac@meta.data <- filtered@meta.data

saveRDS(filtered, paste0(genotype[i], 'filtered.rds'))
saveRDS(filtered.rna, paste0(genotype[i], 'filteredRNA.rds'))
saveRDS(filtered.atac, paste0(genotype[i], 'filteredATAC.rds'))

}


