#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO_DIR/dataWNN_monocle3.rds')
Idents(data) <- 'seurat_clusters'

abT <- subset(data, idents = c('8', '6', '5', '9', '3', '10'), invert = T)
ILC <- subset(data, idents = c('4', '11', '9', '3', '10'))

genes <- c('Rag1', 'Rag2', 'Ptcra', 'Lef1', 'Cd3e', 'Bcl11b', 'Cd3d', 'Cd3g', 'Tcf12')
regions <- c('chr2-101529844-101699843', 'chr2-101529844-101699843', 'chr17-46753051-46771551', 'chr3-130912082-131230082', 
             'chr9-44992417-45018917', 'chr12-106860000-108060000', 'chr9-44967417-45026417', 'chr9-44992417-45018917', 'chr9-71789193-72229193')

for (i in 1:length(genes)){
  
  E2A <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3'])
  E2A$peak <- row.names(E2A)
  E2A <- as.data.frame(E2A[E2A$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3']` == '1', 'peak'])
  colnames(E2A) <- 'peak'
  E2A <- E2A %>% separate(peak, c('chrom', 'start', 'end'), sep='-')
  
  RBPJ <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1'])
  RBPJ$peak <- row.names(RBPJ)
  RBPJ <- as.data.frame(RBPJ[RBPJ$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1']` == '1', 'peak'])
  colnames(RBPJ) <- 'peak'
  RBPJ <- RBPJ %>% separate(peak, c('chrom', 'start', 'end'), sep='-')
  
  p = NULL
  p <- CoveragePlot(subset(data, genotype == 'WT'), region = regions[i], assay = 'ATAC', peaks = F, links = F)
  png(file = paste0(genes[i], '.WT.cover.png'), width = 450, height = 300)
  plot(p)
  dev.off()
  
  region.df <- data.frame(t(data.frame(strsplit(regions[i], '-'))))
  row.names(region.df) = NULL
  colnames(region.df) <- c('chrom', 'start', 'end')
  
  merged.RBPJ <- bedtoolsr::bt.intersect(region.df, RBPJ)
  merged.RBPJ <- paste(merged.RBPJ$V1, merged.RBPJ$V2, merged.RBPJ$V3, sep = '-')
  merged.E2A <- bedtoolsr::bt.intersect(region.df, E2A)
  merged.E2A <- paste(merged.E2A$V1, merged.E2A$V2, merged.E2A$V3, sep = '-')
  
  #Cd3d,g 
  #merged.E2A <- merged.E2A[1:9]
  #merged.RBPJ <- merged.RBPJ[1:6]
  
  abT.TF <- as.data.frame(t(abT@assays[['SCT']][genes[i]]))
  abT.TF$cell <- abT@meta.data$genotype
  abT.TF$pt <- abT@meta.data$monocle
  colnames(abT.TF) <- c('tf', 'cell', 'pt')
  abT.gene <- abT.TF[abT.TF$cell == 'WT', ]
  
  abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.RBPJ]))
  abT.TF$mean <- rowMeans(abT.TF)
  abT.TF$cell <- abT@meta.data$genotype
  abT.TF$pt <- abT@meta.data$monocle
  abT.TF$TF <- 'RBPJ'
  #colnames(abT.TF) <- c('tf', 'cell', 'pt')
  abT.RBPJ <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]
  
  abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.E2A]))
  abT.TF$mean <- rowMeans(abT.TF)
  abT.TF$cell <- abT@meta.data$genotype
  abT.TF$pt <- abT@meta.data$monocle
  abT.TF$TF <- 'E2A'
  #colnames(abT.TF) <- c('tf', 'cell', 'pt')
  abT.E2A <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]
  abT.atac <- rbind(abT.RBPJ, abT.E2A)
  
  ylim.prim <- c(0, mean(abT.gene[abT.gene$pt > 12, 'tf']))
  ylim.sec <- c(0, mean(abT.atac[abT.atac$pt > 12, 'mean']))
  
  
  b <- diff(ylim.prim)/diff(ylim.sec)
  a <- ylim.prim[1] - b*ylim.sec[1]
  
  p = NULL
  p <- ggplot(NULL)
  p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'red3', size = 3)
  p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                       color = 'darkgreen', linetype = 'dashed', size = 3)
  p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                       color = 'blue4', linetype = 'dashed', size = 3)
  p <- p + theme_classic()
  p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ ./b, name = 'ATAC'))
  p
  
  png(file = paste0('RBPJ_E2A.Tcell.', genes[i], '.loci.png'), width = 900, height = 600)
  plot(p)
  dev.off()
  
  
  ILC.TF <- as.data.frame(t(ILC@assays[['SCT']][genes[i]]))
  ILC.TF$cell <- ILC@meta.data$genotype
  ILC.TF$pt <- ILC@meta.data$monocle
  colnames(ILC.TF) <- c('tf', 'cell', 'pt')
  ILC.gene <- ILC.TF[ILC.TF$cell == 'WT', ]
  
  ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.RBPJ]))
  ILC.TF$mean <- rowMeans(ILC.TF)
  ILC.TF$cell <- ILC@meta.data$genotype
  ILC.TF$pt <- ILC@meta.data$monocle
  ILC.TF$TF <- 'RBPJ'
  #colnames(abT.TF) <- c('tf', 'cell', 'pt')
  ILC.RBPJ <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]
  
  ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.E2A]))
  ILC.TF$mean <- rowMeans(ILC.TF)
  ILC.TF$cell <- ILC@meta.data$genotype
  ILC.TF$pt <- ILC@meta.data$monocle
  ILC.TF$TF <- 'E2A'
  #colnames(ILC.TF) <- c('tf', 'cell', 'pt')
  ILC.E2A <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]
  
  ILC.atac <- rbind(ILC.RBPJ, ILC.E2A)
  
  p = NULL
  p <- ggplot(NULL)
  p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'grey20', size = 3)
  p <- p + geom_smooth(data = ILC.gene, aes(pt, tf), color = 'red3', size = 3)
  p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                       color = 'darkgreen', linetype = 'dashed', size = 3)
  p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                       color = 'blue4', linetype = 'dashed', size = 3)
  p <- p + theme_classic()
  p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ ./b, name = 'ATAC'), breaks = seq(0, 4, length = 5))
  p
  
  png(file = paste0('RBPJ_E2A.ILC.', genes[i], '.loci.png'), width = 900, height = 600)
  plot(p)
  dev.off()

}


## Rag2 ##
i = 2

E2A <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3'])
E2A$peak <- row.names(E2A)
E2A <- as.data.frame(E2A[E2A$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA0522.3']` == '1', 'peak'])
colnames(E2A) <- 'peak'
E2A <- E2A %>% separate(peak, c('chrom', 'start', 'end'), sep='-')

RBPJ <- as.data.frame(data@assays[['ATAC']]@motifs@data[,data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1'])
RBPJ$peak <- row.names(RBPJ)
RBPJ <- as.data.frame(RBPJ[RBPJ$`data@assays[['ATAC']]@motifs@data[, data@assays[['ATAC']]@motifs@data@Dimnames[[2]] == 'MA1621.1']` == '1', 'peak'])
colnames(RBPJ) <- 'peak'
RBPJ <- RBPJ %>% separate(peak, c('chrom', 'start', 'end'), sep='-')

region.df <- data.frame(t(data.frame(strsplit(regions[i], '-'))))
row.names(region.df) = NULL
colnames(region.df) <- c('chrom', 'start', 'end')

merged.RBPJ <- bedtoolsr::bt.intersect(region.df, RBPJ)
merged.RBPJ <- paste(merged.RBPJ$V1, merged.RBPJ$V2, merged.RBPJ$V3, sep = '-')
merged.E2A <- bedtoolsr::bt.intersect(region.df, E2A)
merged.E2A <- paste(merged.E2A$V1, merged.E2A$V2, merged.E2A$V3, sep = '-')

abT.TF <- as.data.frame(t(abT@assays[['SCT']][genes[i]]))
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.gene <- abT.TF[abT.TF$cell == 'WT', ]

abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.RBPJ]))
abT.TF$mean <- rowMeans(abT.TF)
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
abT.TF$TF <- 'RBPJ'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.RBPJ <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.E2A]))
abT.TF$mean <- rowMeans(abT.TF)
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
abT.TF$TF <- 'E2A'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.E2A <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]
abT.atac <- rbind(abT.RBPJ, abT.E2A)

ylim.prim <- c(0, mean(abT.gene[abT.gene$pt > 12, 'tf']))
ylim.sec <- c(0, mean(abT.atac[abT.atac$pt > 12, 'mean']))


b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
b <- 2.3

p = NULL
p <- ggplot(NULL)
p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'red3', size = 3)
p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                     color = 'darkgreen', linetype = 'dashed', size = 3)
p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                     color = 'blue4', linetype = 'dashed', size = 3)
p <- p + theme_classic()
p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ ./b, name = 'ATAC'))
p

png(file = paste0('RBPJ_E2A.Tcell.', genes[i], '.loci.v2.png'), width = 900, height = 600)
plot(p)
dev.off()


ILC.TF <- as.data.frame(t(ILC@assays[['SCT']][genes[i]]))
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
colnames(ILC.TF) <- c('tf', 'cell', 'pt')
ILC.gene <- ILC.TF[ILC.TF$cell == 'WT', ]

ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.RBPJ]))
ILC.TF$mean <- rowMeans(ILC.TF)
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
ILC.TF$TF <- 'RBPJ'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
ILC.RBPJ <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.E2A]))
ILC.TF$mean <- rowMeans(ILC.TF)
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
ILC.TF$TF <- 'E2A'
#colnames(ILC.TF) <- c('tf', 'cell', 'pt')
ILC.E2A <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

ILC.atac <- rbind(ILC.RBPJ, ILC.E2A)

p = NULL
p <- ggplot(NULL)
p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'grey20', size = 3)
p <- p + geom_smooth(data = ILC.gene, aes(pt, tf), color = 'red3', size = 3)
p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                     color = 'darkgreen', linetype = 'dashed', size = 3)
p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                     color = 'blue4', linetype = 'dashed', size = 3)
p <- p + theme_classic()
p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ ./b, name = 'ATAC'))

png(file = paste0('RBPJ_E2A.ILC.', genes[i], '.loci.v2.png'), width = 900, height = 600)
plot(p)
dev.off()


## Xccr6 ##
#E2A binding site defined by ChIP (mm10):chr15-82027493-82028448
regions <- 'chr15-81981570-82053570'
genes <- 'Xrcc6'
E2A <- 'chr15-82027493-82028448'
E2A <- data.frame(t(data.frame(strsplit(E2A, '-'))))
row.names(E2A) = NULL
colnames(E2A) <- c('chrom', 'start', 'end')

peaks <- row.names(data@assays[['ATAC']][grep('^chr15', row.names(data@assays[['ATAC']]))])
peaks <- data.frame(t(data.frame(strsplit(peaks, '-'))))
row.names(peaks) = NULL
colnames(peaks) <- c('chrom', 'start', 'end')

E2A <- bedtoolsr::bt.intersect(peaks, E2A)

region.df <- data.frame(t(data.frame(strsplit(regions, '-'))))
row.names(region.df) = NULL
colnames(region.df) <- c('chrom', 'start', 'end')
merged.E2A <- bedtoolsr::bt.intersect(region.df, E2A)
merged.E2A <- paste(merged.E2A$V1, merged.E2A$V2, merged.E2A$V3, sep = '-')

merged.RBPJ <- bedtoolsr::bt.intersect(region.df, RBPJ)
merged.RBPJ <- paste(merged.RBPJ$V1, merged.RBPJ$V2, merged.RBPJ$V3, sep = '-')


abT.TF <- as.data.frame(t(abT@assays[['SCT']][genes]))
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.gene <- abT.TF[abT.TF$cell == 'WT', ]

abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.RBPJ]))
abT.TF$mean <- rowMeans(abT.TF)
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
abT.TF$TF <- 'RBPJ'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.RBPJ <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

abT.TF <- as.data.frame(t(abT@assays[['ATAC']][merged.E2A]))
abT.TF$mean <- rowMeans(abT.TF)
abT.TF$cell <- abT@meta.data$genotype
abT.TF$pt <- abT@meta.data$monocle
abT.TF$TF <- 'E2A'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
abT.E2A <- abT.TF[abT.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

abT.atac <- rbind(abT.RBPJ, abT.E2A)


ylim.prim <- c(0, ceiling(mean(abT.gene[abT.gene$pt > 12, 'tf'])))
ylim.sec <- c(0, round(mean(abT.atac[abT.atac$pt > 12, 'mean']) + 0.05, digits = 2))

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

p = NULL
p <- ggplot(NULL)
p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'red3', size = 3)
p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                     color = 'darkgreen', linetype = 'dashed', size = 3)
p <- p + geom_smooth(data = abT.atac[abT.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                     color = 'blue4', linetype = 'dashed', size = 3)
p <- p + theme_classic()
p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ (. - a)/b, name = 'ATAC'))

png(file = paste0('RBPJ_E2A.Tcell.', genes, '.loci.png'), width = 900, height = 600)
plot(p)
dev.off()


ILC.TF <- as.data.frame(t(ILC@assays[['SCT']][genes]))
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
colnames(ILC.TF) <- c('tf', 'cell', 'pt')
ILC.gene <- ILC.TF[ILC.TF$cell == 'WT', ]

ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.RBPJ]))
ILC.TF$mean <- rowMeans(ILC.TF)
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
ILC.TF$TF <- 'RBPJ'
#colnames(abT.TF) <- c('tf', 'cell', 'pt')
ILC.RBPJ <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

ILC.TF <- as.data.frame(t(ILC@assays[['ATAC']][merged.E2A]))
ILC.TF$mean <- rowMeans(ILC.TF)
ILC.TF$cell <- ILC@meta.data$genotype
ILC.TF$pt <- ILC@meta.data$monocle
ILC.TF$TF <- 'E2A'
#colnames(ILC.TF) <- c('tf', 'cell', 'pt')
ILC.E2A <- ILC.TF[ILC.TF$cell == 'WT', c('mean', 'cell', 'pt', 'TF')]

ILC.atac <- rbind(ILC.RBPJ, ILC.E2A)

p = NULL
p <- ggplot(NULL)
p <- p + geom_smooth(data = abT.gene, aes(pt, tf), color = 'grey20', size = 3)
p <- p + geom_smooth(data = ILC.gene, aes(pt, tf), color = 'red3', size = 3)
p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'RBPJ',], aes(pt, mean*b), 
                     color = 'darkgreen', linetype = 'dashed', size = 3)
p <- p + geom_smooth(data = ILC.atac[ILC.atac$TF %in% 'E2A',], aes(pt, mean*b), 
                     color = 'blue4', linetype = 'dashed', size = 3)
p <- p + theme_classic()
p <- p + scale_y_continuous(limits = c(0, NA), 'RNA', sec.axis = sec_axis(~ (. - a)/b, name = 'ATAC'))

png(file = paste0('RBPJ_E2A.ILC.', genes, '.loci.png'), width = 900, height = 600)
plot(p)
dev.off()


