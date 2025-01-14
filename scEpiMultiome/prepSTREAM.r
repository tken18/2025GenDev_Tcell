#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO_DIR/thymocyteWNN.rds')

# Write an ATAC count matrix -----
data <- subset(data, idents = c('8', '5', '6'), invert = T)
count <- data@assays$ATAC@counts
feature.names.bed <- t(as.data.frame(strsplit(count@Dimnames[[1]], '-')))
row.names(feature.names.bed) = NULL
lympho.bar <- count@Dimnames[[2]]
count@Dimnames[[2]] = NULL
count@Dimnames[[1]] = NULL
writeMM(count, 'atac.lympho.matrix.mtx')
write.table(feature.names.bed, 'atac.lympho.features.bed', sep = '\t', row.names = F, col.names = F, quote = F)
write.table(lympho.bar, 'barcodes.lympho.tsv', sep = '\t', row.names = F, col.names = F, quote = F)

# Write cell labels -----
pal <- c('Set_colors')
cluster.label <- as.array(data@active.ident)
cell.colors <- as.data.frame(sccore::fac2col(cluster.label, level.colors = pal))
cell.colors$cell <- row.names(cell.colors)
cluster.label <- as.data.frame(cluster.label)
cluster.label$cell <- row.names(cluster.label)
write.table(cluster.label, 'label.tsv', col.names = F, row.names = T, sep = '\t', quote = F)
write.table(cell.colors, 'col.tsv', col.names = F, row.names = T, sep = '\t', quote = F)

genotype <- as.data.frame(data@meta.data$genotype)
row.names(genotype) <- row.names(cluster.label)
write.table(genotype, 'genotype.tsv', col.names = F, row.names = T, sep = '\t', quote = F)


# Calculate Z-score of motifss -----
file.count = '/PATH_TO/atac.lympho.matrix.mtx'
file.region = '/PATH_TO/atac.lympho.features.bed'
file.sample = '/PATH_TO/barcodes.lympho.tsv'
dir.output = '/PATH_TO_DIR/'
file.format = 'mtx'
genome = 'BSgenome.Mmusculus.UCSC.mm10'
feature = 'motif'
species = 'Homo sapiens' #Mus musculus
n_jobs = '1'
resize_peak = NULL
peak_width = NULL
resize_peak = FALSE
  
suppressMessages(library(chromVAR))
suppressMessages(library(Matrix))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(BiocParallel))
suppressMessages(library(data.table))
register(MulticoreParam(n_jobs))
peaks <- makeGRangesFromDataFrame(data.frame(fread(file.region,col.names=c('seqnames','start','end'))))  
if(resize_peak==FALSE){
  print('none')
}else{
  print('resize regions/peaks ...')
  peaks <- resize(peaks, width = peak_width, fix = 'center')
}
samples <- data.frame(fread(file.sample,header=FALSE))
# Import counts
if(file.format=='mtx'){
  counts = readMM(file.count)
  colnames(counts) <- samples[,1]
}else{
  m <- data.matrix(data.frame(fread(cmd=paste0('zcat < ', file.count))))
  counts <- sparseMatrix(i = m[,1], j = m[,2], x=m[,3])
  colnames(counts) <- samples[,1]
}
# Make RangedSummarizedExperiment
SE <- SummarizedExperiment(
  rowRanges = peaks,
  colData = samples,
  assays = list(counts = counts)
)
# Add GC bias
SE <- addGCBias(SE, genome = genome)
SE <- filterPeaks(SE, non_overlapping = TRUE)
bg <- getBackgroundPeaks(SE)
# compute motif deviations
suppressMessages(library('JASPAR2020')) 
suppressMessages(library(motifmatchr))
motifs = getJasparMotifs(species = species)    
motif_ix <- matchMotifs(motifs, SE, genome = genome)
dev <- computeDeviations(object = SE, annotations = motif_ix, background_peaks = bg)  
devTable <- assays(dev)[['deviations']]
write.table(devTable, file = gzfile(file.path(dir.output,'zscores.tsv.gz')), sep = '\t', quote = FALSE, row.names = TRUE, col.names=NA)

