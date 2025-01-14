#!/bin/Rscript

setwd('/WD_PATH/')


# Read data -----
data <- readRDS('/PATH_TO_DIR/thymocyteWNN_monocle3.rds')

Idents(data) <- 'seurat_clusters'
abT <- subset(data, idents = c('8', '6', '5', '9', '3', '10'), invert = T)
ILC <- subset(data, idents = c('4', '11', '9', '3', '10'))


# Make a motif list -----
motif.list <- read.table('/PATH_TO_DIR/JASPAR2020_CORE_redundant_pfms_jaspar.txt', header = F, sep = '-')
motif.list <- sub('>', '', motif.list[grep('^>MA', motif.list$V1),])
motif.list <- as.data.frame(t(as.data.frame(str_split(motif.list, '\t'))))
motif.list <- motif.list[motif.list$V1 %in% intersect(motif.list$V1, row.names(data@assays$chromvar)),]
TFs <- motif.list[grep('Rbpj', motif.list$V2, ignore.case = T),]
for (gene in TF.list[2:length(TF.list)]){
  
  tmp <- motif.list[grep(gene, motif.list$V2, ignore.case = T),]
  TFs <- rbind(TFs, tmp)
  
}
TFs.motif <- TFs$V1


# Wave plots of motif activity -----
for (i in 1:length(TFs.motif)){
  
  abT.TF <- as.data.frame(t(abT@assays[['chromvar']][TFs.motif[i]]))
  abT.TF$cell <- abT@meta.data$genotype
  abT.TF$pt <- abT@meta.data$monocle
  colnames(abT.TF) <- c('tf.motif', 'cell', 'pt')
  abT.WT <- abT.TF[abT.TF$cell == 'WT', ]
  abT.KO <- abT.TF[abT.TF$cell == 'dKO', ]
  
  ILC.TF <- as.data.frame(t(ILC@assays[['chromvar']][TFs.motif[i]]))
  ILC.TF$cell <- ILC@meta.data$genotype
  ILC.TF$pt <- ILC@meta.data$monocle
  colnames(ILC.TF) <- c('tf.motif', 'cell', 'pt')
  ILC.WT <- ILC.TF[ILC.TF$cell == 'WT', ]
  ILC.KO <- ILC.TF[ILC.TF$cell == 'dKO', ]
  
  p = NULL
  p <- ggplot(NULL)
  p <- p + geom_smooth(data = abT.WT, aes(pt, tf.motif), color = 'blue3', size = 3)
  p <- p + geom_smooth(data = abT.KO, aes(pt, tf.motif), color = 'chocolate2', size = 3)
  p <- p + geom_smooth(data = ILC.WT, aes(pt, tf.motif), color = 'blue3', linetype = 'dashed', size = 3)
  p <- p + geom_smooth(data = ILC.KO, aes(pt, tf.motif), color = 'chocolate2', linetype = 'dashed', size = 3)
  p <- p + theme_classic()
  p
  
  pdf(file=paste0('/PATH_TO/motif_pdf/', TFs.motif[i], '_', TFs[TFs$V1 == TFs.motif[i], 'V2'], '.pdf'), width = 9, height = 6)
  plot(p)
  dev.off()
  
  png(file=paste0('/PATH_TO/motif_png/', TFs.motif[i], '_', TFs[TFs$V1 == TFs.motif[i], 'V2'], '.png'), width = 900, height = 600)
  plot(p)
  dev.off()
  
}


# Wave plots of TF expression -----
TF.list <- read.table('TF.list.txt', header = F, sep = '\t')
TF.list <- intersect(TF.list$V1, row.names(data@assays$RNA))

for (i in TF.list){
  
  abT.TF <- as.data.frame(t(abT@assays[['SCT']][i]))
  abT.TF$cell <- abT@meta.data$genotype
  abT.TF$pt <- abT@meta.data$monocle
  colnames(abT.TF) <- c('tf.gene', 'cell', 'pt')
  abT.WT <- abT.TF[abT.TF$cell == 'WT', ]
  abT.KO <- abT.TF[abT.TF$cell == 'dKO', ]
  
  ILC.TF <- as.data.frame(t(ILC@assays[['SCT']][i]))
  ILC.TF$cell <- ILC@meta.data$genotype
  ILC.TF$pt <- ILC@meta.data$monocle
  colnames(ILC.TF) <- c('tf.gene', 'cell', 'pt')
  ILC.WT <- ILC.TF[ILC.TF$cell == 'WT', ]
  ILC.KO <- ILC.TF[ILC.TF$cell == 'dKO', ]
  
  p = NULL
  p <- ggplot(NULL)
  p <- p + geom_smooth(data = abT.WT, aes(pt, tf.gene), color = 'blue3', size = 3)
  p <- p + geom_smooth(data = abT.KO, aes(pt, tf.gene), color = 'chocolate2', size = 3)
  p <- p + geom_smooth(data = ILC.WT, aes(pt, tf.gene), color = 'blue3', linetype = 'dashed', size = 3)
  p <- p + geom_smooth(data = ILC.KO, aes(pt, tf.gene), color = 'chocolate2', linetype = 'dashed', size = 3)
  p <- p + theme_classic()
  p
  
  pdf(file=paste0('/PATH_TO/gene_pdf/', i, '.pdf'), width = 9, height = 6)
  plot(p)
  dev.off()
  
  png(file=paste0('/PATH_TO/gene_png/', i, '.png'), width = 900, height = 600)
  plot(p)
  dev.off()
  
}


# Density plots -----
abT.TF$cluster <- abT@active.ident
ILC.TF$cluster <- ILC@active.ident

library(plyr)
abT.TF.mu <- ddply(abT.TF, 'cluster', summarise, grp.mean = mean(pt))
ggplot(abT.TF, aes(x = pt, color = cluster)) +
  geom_density(size = 2) + theme_classic() + scale_color_manual(values=c('#EE0011FF', '#0C5BB0FF', '#15983DFF', '#FA6B09FF', '#FEC10BFF', '#6351A0FF'))

ILC.TF.mu <- ddply(ILC.TF, 'cluster', summarise, grp.mean = mean(pt))
ggplot(ILC.TF, aes(x = pt, color = cluster)) +
  geom_density(size = 2) + theme_classic() + scale_color_manual(values=c('#EC579AFF', '#FA6B09FF', '#9A703EFF', '#F5BACFFF', '#6351A0FF'))



