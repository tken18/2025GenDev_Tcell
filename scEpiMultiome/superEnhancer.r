#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO_DATA/integrated.data.rds')

# Read core super enhancer -----
Tcore.bed <- read.table('/PATH_TO_BED/T.bed', sep = '\t', header = F)
Tcore.bed <- paste(Tcore.bed$V4, Tcore.bed$V5, Tcore.bed$V6, sep = '-')
ILCcore.bed <- read.table('/PATH_TO_BED/ILC2.bed', sep = '\t', header = F)
ILCcore.bed <- paste(ILCcore.bed$V4, ILCcore.bed$V5, ILCcore.bed$V6, sep = '-')
peaks <- data@assays$ATAC@data
DefaultAssay(data) <- 'ATAC'


# Analysis core super enhancer -----
Tcore.peaks <- as.data.frame(t(peaks[peaks@Dimnames[[1]] %in% Tcore.bed,]))
Tcore.peaks$means <- rowMeans(Tcore.peaks)
Tcore.peaks$cluster <- data@meta.data$seurat_clusters
Tcore.peaks$genotype <- data@meta.data$genotype
Tcore.peaks$cluster <- factor(Tcore.peaks$cluster, levels = c('4', '11', '2', '1', '7', '0', '9', '3', '10', '5', '6', '8'))
Tcore.peaks$group <- paste(Tcore.peaks$cluster, Tcore.peaks$genotype, sep = '.')

ILCcore.peaks <- as.data.frame(t(peaks[peaks@Dimnames[[1]] %in% ILCcore.bed,]))
ILCcore.peaks$means <- rowMeans(ILCcore.peaks)
ILCcore.peaks$cluster <- data@meta.data$seurat_clusters
ILCcore.peaks$genotype <- data@meta.data$genotype
ILCcore.peaks$cluster <- factor(ILCcore.peaks$cluster, levels = c('4', '11', '2', '1', '7', '0', '9', '3', '10', '5', '6', '8'))
ILCcore.peaks$group <- paste(ILCcore.peaks$cluster, ILCcore.peaks$genotype, sep = '.')

ggplot() + theme_classic() +
  geom_point(data = Tcore.peaks[Tcore.peaks$cluster != '8' & Tcore.peaks$cluster != '5' & Tcore.peaks$cluster != '6',],
             aes(x = cluster, y = means, color = genotype, fill = cluster),
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 0.6) +
  geom_boxplot(data = Tcore.peaks[Tcore.peaks$cluster != '8' & Tcore.peaks$cluster != '5' & Tcore.peaks$cluster != '6',],
               aes(x = cluster, y = means, color = genotype, fill = cluster),
               outlier.colour = NA, position = position_dodge(width = 0.8), width=0.6) +
  scale_fill_manual(values = pallet_ad_white) +
  scale_color_manual(values = pallet_ad_man)


ggplot() + theme_classic() +
  geom_point(data = ILCcore.peaks[ILCcore.peaks$cluster != '8' & ILCcore.peaks$cluster != '5' & ILCcore.peaks$cluster != '6',],
             aes(x = cluster, y = means, color = genotype, fill = cluster),
             position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 0.6) +
  geom_boxplot(data = ILCcore.peaks[ILCcore.peaks$cluster != '8' & ILCcore.peaks$cluster != '5' & ILCcore.peaks$cluster != '6',],
               aes(x = cluster, y = means, color = genotype, fill = cluster),
               outlier.colour = NA, position = position_dodge(width = 0.8), width=0.6) +
  scale_fill_manual(values = pallet_ad_white) +
  scale_color_manual(values = pallet_ad_man)


res_Table = NULL
for (i in as.character(0:11)){
  
  WT <- Tcore.peaks[Tcore.peaks$group %in% paste0(i, '.WT'), 'means']
  dKO <- Tcore.peaks[Tcore.peaks$group %in% paste0(i, '.dKO'), 'means']
  T.core.test <- wilcox.test(WT, dKO)
  t <- T.core.test$statistic
  p <- T.core.test$p.value
  tmp <- data.frame(i, t, p)
  res_Table <- rbind(res_Table, tmp)
  
}
write.table(res_Table, 'TcoreTestPerCell.txt', quote = F, sep = '\t', row.names = F)


res_Table = NULL
for (i in as.character(0:11)){
  
  WT <- ILCcore.peaks[ILCcore.peaks$group %in% paste0(i, '.WT'), 'means']
  dKO <- ILCcore.peaks[ILCcore.peaks$group %in% paste0(i, '.dKO'), 'means']
  ILC.core.test <- wilcox.test(WT, dKO)
  t <- ILC.core.test$statistic
  p <- ILC.core.test$p.value
  tmp <- data.frame(i, t, p)
  res_Table <- rbind(res_Table, tmp)
  
}
write.table(res_Table, 'ILC2coreTestPerCell.txt', quote = F, sep = '\t', row.names = F)


