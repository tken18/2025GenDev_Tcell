#!/bin/Rscript

setwd('/WD_PATH/')

data <- readRDS('/PATH_TO/thymocyteWNN_monocle3.rds')
transient.cds <- readRDS('/PATH_TO/subsetTransientCells.rds')

Tcell <- subset(data, idents = '8', genotype == 'WT')
ILCp <- subset(data, idents = '9', genotype == 'WT')
transient <- subset(data, cells = colnames(transient.cds), genotype == 'dKO')


TFs <- c('Tcf3', 'Gata3', 'Tcf7', 'Runx1', 'Runx2', 'Runx3', 'Rora', 'Notch1', 'Spi1', 'Nfil3', 'Erg', 'Ets1')
for (i in TFs){
  abT.TF <- as.data.frame(t(Tcell@assays[["SCT"]][i]))
  abT.TF$cell <- 'Tcell'
  colnames(abT.TF) <- c('tf', 'cell')
  
  ILC.TF <- as.data.frame(t(ILCp@assays[["SCT"]][i]))
  ILC.TF$cell <- 'ILCp'
  colnames(ILC.TF) <- c('tf', 'cell')
  
  transient.TF <- as.data.frame(t(transient@assays[["SCT"]][i]))
  transient.TF$cell <- 'transient'
  colnames(transient.TF) <- c('tf', 'cell')
  
  tmp <- rbind(abT.TF, ILC.TF, transient.TF)
  tmp$cell <- factor(tmp$cell, levels = c('Tcell', 'transient', 'ILCp'))
  
  test <- pairwise.t.test(tmp$tf, tmp$cell, p.adj = "bonf")
  write.table(test$p.value, paste0('test/', i, '.t.test.txt'), quote = F)
  
  p <- ggplot(NULL) + theme_classic()
  
  p <- p + geom_point(data = tmp,
               aes(x = cell, y = tf, fill = cell),
               position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 0.6)
  
  p <- p + geom_boxplot(data = tmp,
                 aes(x = cell, y = tf, fill = cell),
                 outlier.colour = NA, position = position_dodge(width = 0.8), width=0.6)
  
  png(paste0(i, '.png'), height = 500, width = 600)
  plot(p)
  dev.off()
  
}



TFs <- c('MA0522.3', 'MA0037.3', 'MA0769.2', 'MA0002.2', 'MA0511.2', 'MA0684.2', 
         'MA0071.1', 'MA1116.1', 'MA0080.5', 'MA0025.2', 'MA0098.3', 'MA0474.2')
for (i in TFs){
  
  abT.TF <- as.data.frame(t(Tcell@assays[["chromvar"]][i]))
  abT.TF$cell <- 'Tcell'
  colnames(abT.TF) <- c('tf', 'cell')
  
  ILC.TF <- as.data.frame(t(ILCp@assays[["chromvar"]][i]))
  ILC.TF$cell <- 'ILCp'
  colnames(ILC.TF) <- c('tf', 'cell')
  
  transient.TF <- as.data.frame(t(transient@assays[["chromvar"]][i]))
  transient.TF$cell <- 'transient'
  colnames(transient.TF) <- c('tf', 'cell')
  
  tmp <- rbind(abT.TF, ILC.TF, transient.TF)
  tmp$cell <- factor(tmp$cell, levels = c('Tcell', 'transient', 'ILCp'))
  
  test <- pairwise.t.test(tmp$tf, tmp$cell, p.adj = "bonf")
  write.table(test$p.value, paste0('test/', i, '.t.test.txt'), quote = F)
  
  p <- ggplot(NULL) + theme_classic()
  
  p <- p + geom_point(data = tmp,
                      aes(x = cell, y = tf, fill = cell),
                      position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.8), size = 0.6)
  
  p <- p + geom_boxplot(data = tmp,
                 aes(x = cell, y = tf, fill = cell),
                 outlier.colour = NA, position = position_dodge(width = 0.8), width=0.6)
  
  png(paste0(i, '.png'), height = 500, width = 600)
  plot(p)
  dev.off()
  
}