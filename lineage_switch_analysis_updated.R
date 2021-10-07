source('scDataAnalysis_Utilities.R')
plotExprFrac <- function(seurat.obj, feature = 'MPO', assay = 'RNA',
                         group.name = 'sample',
                         plot.group.order){
  feature.exprs = seurat.obj[[assay]]@counts[feature, ]
  g.types = unique(seurat.obj@meta.data[, group.name])
  exprs.g <- lapply(g.types, function(x) feature.exprs[seurat.obj[[group.name]] == x] )
  names(exprs.g) = g.types
  
  dd.dt = data.table('group' = g.types,
                     'frac' = sapply(g.types, function(x) mean(exprs.g[[x]] > 0)))
  dd.dt$group = factor(dd.dt$group, levels = plot.group.order)
  dd.dt$frac = round(dd.dt$frac, 3)
  ymax = min(1, max(dd.dt$frac)+0.05)
  p0 <- ggplot(data = dd.dt, aes(x = group, y = frac, fill = group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab('Fraction of expression') + ylim(0, ymax) +
    geom_text(aes(label = frac, x = group, y = frac), 
              position = position_dodge(width = 0.8), vjust = -0.6) 
  return(p0)
  
}

plotExprFrac1 <- function(seurat.obj, feature = 'MPO', assay = 'RNA',
                         group.name = 'sample',
                         plot.group.order,
                         do.stat.test = FALSE,
                         show.pv = FALSE,
                         test.control.group = 'HD'){
  feature.exprs = seurat.obj[[assay]]@counts[feature, ]
  g.types = unique(seurat.obj@meta.data[, group.name])
  exprs.g <- lapply(g.types, function(x) feature.exprs[seurat.obj[[group.name]] == x] )
  names(exprs.g) = g.types
  pvs = rep(1, length(g.types))
  names(pvs) = g.types
  if(do.stat.test){
    for(type0 in g.types){
      if(type0 == test.control.group) next
      pvs[type0] = pbinom(sum(exprs.g[[type0]] > 0), 
                          length(exprs.g[[type0]]),
                          prob = mean(exprs.g[[test.control.group]] > 0),
                          lower.tail = F)
      pvs[type0] = min(pvs[type0], 1-pvs[type0])
    }
  }
  
  
  if(show.pv){
    dd.dt = data.table('group' = g.types,
                       'pv' = paste0('p=', format(pvs, digit = 2)),
                       'frac' = sapply(g.types, function(x) mean(exprs.g[[x]] > 0)))
    dd.dt[group == test.control.group]$pv = 'background'
    dd.dt$group = factor(dd.dt$group, levels = plot.group.order)
    dd.dt$frac = round(dd.dt$frac, 3)
    ymax = min(1, max(dd.dt$frac) * 1.1)
    
    p0 <- ggplot(data = dd.dt, aes(x = group, y = frac, fill = group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      ylab('Fraction of expression') + ylim(0, ymax) +
      geom_text(aes(label = pv, x = group, y = frac), 
                position = position_dodge(width = 0.8), vjust = -0.6) 
    
  }else{
    dd.dt = data.table('group' = g.types,
                       'frac' = sapply(g.types, function(x) mean(exprs.g[[x]] > 0)))
    dd.dt$group = factor(dd.dt$group, levels = plot.group.order)
    dd.dt$frac = round(dd.dt$frac, 3)
    ymax = min(1, max(dd.dt$frac) * 1.1)
    p0 <- ggplot(data = dd.dt, aes(x = group, y = frac, fill = group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      ylab('Fraction of expression') + ylim(0, ymax) +
      geom_text(aes(label = frac, x = group, y = frac), 
                position = position_dodge(width = 0.8), vjust = -0.6) 
    
    
  }
  return(list('plot' = p0, 'pvs' = pvs))
  
}


## plot myeloid lineage percentage -- updated and with LMPP ####
seurat.rna <- readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.nkt <- subset(seurat.rna, Ctype0 %in% c('T', 'NKT') & Ctype_updated != 'doublets')

seurat.blast = subset(seurat.rna, Ctype0 == 'Blasts')
rm(seurat.rna)

seurat.atac = readRDS('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap10000.rds')
seurat.nkt.atac = subset(seurat.atac, Ctype0 == 'T/NK')

seurat.blast.atac = subset(seurat.atac, Ctype0 == 'Blasts')
rm(seurat.atac)

getPalette = colorRampPalette(brewer.pal(9, "Paired"))
color_sample = getPalette(18)
names(color_sample) = c( "PAYWJZ", "PAZGKI", "PAYUZM", "1154",  
                         "PAYKGI", "PAZBSZ", "PAYWKL", "PAYSBA", 
                         "PAZBLA", "PAYZLC", "PAYLNH", "PAYUZJ", 
                         "PAYYBG", "PAZFPH", "PAYYNY", "PAZBGV",
                         "PAYZVY", "PAYZWN ")

b.types = c('CLP', 'Pre-pro-B', 'Pro-B',
            'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B')
b.types1 = paste0(b.types, '-like')

m.types = c('pDC', 'cDC', 'DC-Progenitor', 
            'GMP', 'Mono')
m.types1 = paste0(m.types, '-like')
seurat.blast$lineage = 'other'
seurat.blast$lineage <- ifelse(seurat.blast$Ctype_updated %in% m.types1,
                               'M-lineage', seurat.blast$lineage)
seurat.blast$lineage <- ifelse(seurat.blast$Ctype_updated %in% b.types1, 
                               'B-lineage', seurat.blast$lineage)
seurat.blast$lineage <- ifelse(seurat.blast$Ctype_updated == 'LMPP-like', 
                               'LMPP', seurat.blast$lineage)


seurat.blast.atac$lineage = 'other'
seurat.blast.atac$lineage <- ifelse(seurat.blast.atac$projCtype %in% m.types1,
                                    'M-lineage', seurat.blast.atac$lineage)
seurat.blast.atac$lineage <- ifelse(seurat.blast.atac$projCtype %in% b.types1, 
                                    'B-lineage', seurat.blast.atac$lineage)
seurat.blast.atac$lineage <- ifelse(seurat.blast.atac$projCtype == 'LMPP-like', 
                                    'LMPP', seurat.blast.atac$lineage)


## < plot m lineage and LMPP by sample for rna ####
ldata <- data.table('sample' = seurat.blast$sample,
                    'lineage' = seurat.blast$lineage)
ldata.freq <- lapply(names(color_sample), function(x) table(ldata[sample ==x]$lineage)/nrow(ldata[sample ==x]) )
names(ldata.freq) = names(color_sample)
m.freq <- sapply(names(color_sample), function(x) ldata.freq[[x]]['M-lineage'])
lmpp.freq <- sapply(names(color_sample), function(x) ldata.freq[[x]]['LMPP'])

m.freq[is.na(m.freq)] = 0
lmpp.freq[is.na(lmpp.freq)] = 0

pdata <- data.table('sample' = names(color_sample),
                    'freq.m' =  m.freq,
                    'freq.lmpp' = lmpp.freq,
                    'lineage' = rep('Myeloid', length(color_sample)))
pdata = pdata[order(-freq.m)]
pdata$sample = factor(pdata$sample, levels = pdata$sample)

p1 <- ggplot(data = pdata, aes(x = sample, y = freq.m, fill= lineage)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('#1F78B4')) +
  ylab('Fraction of myeloid blasts') + xlab('') + 
  geom_hline(yintercept = 0.02, linetype = 2) +
  theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  annotate('text', x = 16, y = 0.03, label = '2%')
ggsave(p1, filename = 'Figures/summary/LineageSwitch/frac_myeloid_blast_rna_updated.eps', 
       device = 'eps', width = 10, height = 4)

p0 <- ggplot(data = pdata, aes(x = sample, y = freq.lmpp, fill= lineage)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('#33A02C')) +
  ylab('Fraction of myeloid blasts') + xlab('') + 
  geom_hline(yintercept = 0.02, linetype = 2) +
  theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  annotate('text', x = 16, y = 0.03, label = '2%')
ggsave(p0, filename = 'Figures/summary/LineageSwitch/frac_LMPP_blast_rna_updated.eps', 
       device = 'eps', width = 10, height = 4)

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_rna.eps',
           width = 6, height = 6)
cr = cor(pdata$freq.lmpp, pdata$freq.m, method = 'spearman')
cor.test(pdata$freq.lmpp, pdata$freq.m,  method = "spearman")
plot(pdata$freq.lmpp, pdata$freq.m, xlab = 'frac of LMPP', 
     ylab = 'frac of M-lineage', pch = 19)
abline(a = 0, b = 1)
text(x = 0.06, y = 0.3, labels = paste0('r = ', round(cr, 3)))
dev.off()


saveRDS(pdata[, 1:3], file = 'MetaData/M_lineage_LMPP_frac_18MLLr_rna.rds')

## < plot m lineage by sample for atac ####
ldata <- data.table('sample' = seurat.blast.atac$sample,
                    'lineage' = seurat.blast.atac$lineage)
ldata.freq <- lapply(names(color_sample), function(x) table(ldata[sample ==x]$lineage)/nrow(ldata[sample ==x]) )
names(ldata.freq) = names(color_sample)
m.freq <- sapply(names(color_sample), function(x) ldata.freq[[x]]['M-lineage'])
lmpp.freq <- sapply(names(color_sample), function(x) ldata.freq[[x]]['LMPP'])

m.freq[is.na(m.freq)] = 0
lmpp.freq[is.na(lmpp.freq)] = 0

pdata <- data.table('sample' = names(color_sample),
                    'freq.m' =  m.freq,
                    'freq.lmpp' = lmpp.freq,
                    'lineage' = rep('Myeloid', length(color_sample)))
pdata = pdata[order(-freq.m)]
pdata$sample = factor(pdata$sample, levels = pdata$sample)

#pdata$sample = factor(pdata$sample, levels = names(color_sample))

p2 <- ggplot(data = pdata, aes(x = sample, y = freq.m, fill= lineage)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('#1F78B4')) +
  ylab('Fraction of myeloid blasts') + xlab('') + 
  geom_hline(yintercept = 0.02, linetype = 2) +
  theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  annotate('text', x = 16, y = 0.04, label = '2%')
ggsave(p2, filename = 'Figures/summary/LineageSwitch/frac_myeloid_blast_atac_updated.eps', 
       device = 'eps', width = 10, height = 4)

p3 <- ggplot(data = pdata, aes(x = sample, y = freq.lmpp, fill= lineage)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('#33A02C')) +
  ylab('Fraction of myeloid blasts') + xlab('') + 
  theme_classic() + NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) 
ggsave(p3, filename = 'Figures/summary/LineageSwitch/frac_LMPP_blast_atac_updated.eps', 
       device = 'eps', width = 10, height = 4)

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_atac.eps',
           width = 6, height = 6)
cr = cor(pdata$freq.lmpp, pdata$freq.m, method = 'spearman')
cor.test(pdata$freq.lmpp, pdata$freq.m, method = "spearman")
plot(pdata$freq.lmpp, pdata$freq.m, xlab = 'frac of LMPP', 
     ylab = 'frac of M-lineage', pch = 19)
abline(a = 0, b = 1)
text(x = 0.003, y = 0.6, labels = paste0('r = ', round(cr, 3)))
dev.off()

saveRDS(pdata[, 1:3], file = 'MetaData/M_lineage_LMPP_frac_18MLLr_atac.rds')


## plot correlation with the two lineage switch examples ####
## < rna ####
cell.map.rna = readRDS(file = 'MetaData/M_lineage_LMPP_frac_lineageSwitch_rna.rds')
pdata.rna = readRDS('MetaData/M_lineage_LMPP_frac_18MLLr_rna.rds')
pdata0.rna = cell.map.rna[ctype_patient %in% c('HTAN2184_Blasts', 'HTAN1979_Blasts')]
pdata0.rna = subset(pdata0.rna, select = c(ctype_patient, lineage, frac))
pdata0.rna[, 'sample' := gsub('_Blasts', '', ctype_patient), by = ctype_patient]
pdata.rna = rbind(pdata.rna, data.table(sample = 'HTAN2184', freq.m = 0.1561,
                                        freq.lmpp = 0.0183))
pdata.rna = rbind(pdata.rna, data.table(sample = 'HTAN1979', freq.m = 0.0354,
                                        freq.lmpp = 0.0545))

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_rna.eps',
           width = 6, height = 6)
cr = cor(pdata.rna$freq.lmpp, pdata.rna$freq.m, method = 'spearman')
cor.test(pdata.rna$freq.lmpp, pdata.rna$freq.m, method = 'spearman')
plot(pdata.rna$freq.lmpp, pdata.rna$freq.m, xlab = 'frac of LMPP', 
     ylab = 'frac of M-lineage', pch = 19)
lines(x = pdata.rna$freq.lmpp, 
      y = lm(pdata.rna$freq.m ~ pdata.rna$freq.lmpp)$fit)
text(x = 0.06, y = 0.3, labels = paste0('r = ', round(cr, 3)))
dev.off()

pdata.rna$freq.lmpp = -log10(pdata.rna$freq.lmpp + 0.00001)
pdata.rna$freq.m = -log10(pdata.rna$freq.m + 0.00001)

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_rna_log10.eps',
           width = 6, height = 6)
plot(pdata.rna$freq.lmpp, pdata.rna$freq.m, 
     ylab = 'frac of M-lineage(-log10 scale)', 
     pch = 19, xlab = "frac of LMPP(-log10 scale)")
lines(x = pdata.rna$freq.lmpp, 
      y = lm(pdata.rna$freq.m ~ pdata.rna$freq.lmpp)$fit)
text(x = 1.5, y = 4, labels = paste0('r = ', round(cr, 3)))
dev.off()



## < atac ####
cell.map.atac = readRDS(file = 'MetaData/M_lineage_LMPP_frac_lineageSwitch_atac.rds')
pdata.atac = readRDS('MetaData/M_lineage_LMPP_frac_18MLLr_atac.rds')
pdata0.atac = cell.map.atac[ctype_patient %in% c('HTAN2184_Blasts', 'HTAN1979_Blasts')]
pdata0.atac = subset(pdata0.atac, select = c(ctype_patient, lineage, frac))
pdata0.atac[, 'sample' := gsub('_Blasts', '', ctype_patient), by = ctype_patient]
pdata.atac = rbind(pdata.atac, data.table(sample = 'HTAN2184', freq.m = 0.1561,
                                        freq.lmpp = 0.0183))
pdata.atac = rbind(pdata.atac, data.table(sample = 'HTAN1979', freq.m = 0.0354,
                                        freq.lmpp = 0.0545))

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_atac.eps',
           width = 6, height = 6)
cr = cor(pdata.atac$freq.lmpp, pdata.atac$freq.m, method = 'spearman')
cor.test(pdata.atac$freq.lmpp, pdata.atac$freq.m, method = 'spearman')
plot(pdata.atac$freq.lmpp, pdata.atac$freq.m, xlab = 'frac of LMPP', 
     ylab = 'frac of M-lineage', pch = 19)
lines(x = pdata.atac$freq.lmpp, 
      y = lm(pdata.atac$freq.m ~ pdata.atac$freq.lmpp)$fit)
text(x = 0.04, y = 0.3, labels = paste0('r = ', round(cr, 3)))
dev.off()

pdata.atac$freq.lmpp = -log10(pdata.atac$freq.lmpp + 0.00001)
pdata.atac$freq.m = -log10(pdata.atac$freq.m+0.00001)

postscript(file = 'Figures/summary/LineageSwitch/scatter_plot_mlineage_lmpp_frac_atac_log10.eps',
           width = 6, height = 6)
plot(pdata.atac$freq.lmpp, pdata.atac$freq.m, 
     ylab = 'frac of M-lineage(-log10 scale)', 
     pch = 19, xlab = "frac of LMPP(-log10 scale)")
lines(x = pdata.atac$freq.lmpp, 
      y = lm(pdata.atac$freq.m ~ pdata.atac$freq.lmpp)$fit)
text(x = 2, y = 2.5, labels = paste0('r = ', round(cr, 3)))
dev.off()


