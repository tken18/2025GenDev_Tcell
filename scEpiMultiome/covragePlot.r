#!/bin/Rscript

setwd('/WD_PATH/')


# Read data -----
data <- readRDS('/PATH_TO/dataWNN_monocle3.rds')
WT <- readRDS('/PATH_TO/dataWNN_monocle3WT.rds')
WT@active.ident <- factor(WT@active.ident, levels = c('4', '11', '2', '7', '0', '1', '9', '3', '10', '5', '6', '8'))
Idents(data) <- 'seurat_clusters'


# Make motif region lists -----
## E2A ##
E2A <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3'])
E2A$peak <- row.names(E2A)
E2A <- as.data.frame(E2A[E2A$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3']` == '1', 'peak'])
colnames(E2A) <- 'peak'
E2A <- E2A %>% separate(peak, c('chrom', 'start', 'end'), sep = '-')

## RBPJ ##
RBPJ <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1'])
RBPJ$peak <- row.names(RBPJ)
RBPJ <- as.data.frame(RBPJ[RBPJ$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1']` == '1', 'peak'])
colnames(RBPJ) <- 'peak'
RBPJ <- RBPJ %>% separate(peak, c('chrom', 'start', 'end'), sep='-')


# Coverage plots -----
regions <- c('chr2-101529844-101699843', 'chr2-101529844-101699843', 'chr15-81981570-82053570', 'chr17-46753051-46771551', 'chr17-32107055-32176055', 
             'chr2-26446480-26568480', 'chr3-130912082-131230082', 'chr9-44967417-45026417', 'chr9-44967417-45026417', 
             'chr9-44992417-45018917', 'chr9-71789193-72229193', 'chr12-106860000-108060000')
genes <- c('Rag1', 'Rag2', 'Xrcc6', 'Ptcra', 'Notch3', 'Notch1', 'Lef1', 'Cd3g', 'Cd3d', 'Cd3e', 'Tcf12', 'Bcl11b')
yaxis <- c(130, 130, 270, 230, 160, 160, 160, 220, 220, 160, 220, 130)

## E2A ##
for (i in 1:length(genes)){

  region.df <- data.frame(t(data.frame(strsplit(regions[i], '-'))))
  row.names(region.df) = NULL
  colnames(region.df) <- c('chrom', 'start', 'end')
  merged <- bedtoolsr::bt.intersect(region.df, E2A)
  merged <- paste(merged$V1, merged$V2, merged$V3, sep = '-')
  #Cd3d,g merge[1:9]
  #merged <- merged[1:9]
  range.show = StringToGRanges(regions = merged)
  
  p = NULL
  p <- CoveragePlot(WT, region = regions[i], assay = 'ATAC', peaks = F, ranges = range.show, ymax = yaxis[i], 
             region.highlight = range.show, ranges.title = 'E2A', idents = c('4', '11', '2', '7', '0', '1', '9', '3', '10'), links = F) & 
             scale_fill_manual(values = pallet_ad_white)
  png(file = paste0(genes[i], '.WT.cover.E2A.png'), width = 450, height = 300)
  plot(p)
  dev.off()

}

## RBPJ ##
for (i in 1:length(genes)){
  
  region.df <- data.frame(t(data.frame(strsplit(regions[i], '-'))))
  row.names(region.df) = NULL
  colnames(region.df) <- c('chrom', 'start', 'end')
  merged <- bedtoolsr::bt.intersect(region.df, RBPJ)
  merged <- paste(merged$V1, merged$V2, merged$V3, sep = '-')
  #Cd3d,g merge[1:6]
  #merged <- merged[1:6]
  range.show = StringToGRanges(regions = merged)
  
  p = NULL
  p <- CoveragePlot(WT, region = regions[i], assay = 'ATAC', peaks = F, ranges = range.show, ymax = yaxis[i], 
                    region.highlight = range.show, ranges.title = 'RBPJ', idents = c('4', '11', '2', '7', '0', '1', '9', '3', '10'), links = F) & 
    scale_fill_manual(values = pallet_ad_white)
  png(file = paste0(genes[i], '.WT.cover.RBPJ.png'), width = 450, height = 300)
  plot(p)
  dev.off()

}
