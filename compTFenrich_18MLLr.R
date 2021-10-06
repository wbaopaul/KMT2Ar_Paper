## compare TF enrichment between patients and HD

source('scDataAnalysis_Utilities.R')


filterTFenrichByExpr <- function(tf.diff.res, expr.mtx,
                                 frac.thr = 0.1){
  cls = unique(tf.diff.res$cluster1)
  res = NULL
  for(cl0 in cls){
    tf.diff.res0 = tf.diff.res[cluster1 == cl0]
    frac0 = expr.mtx[, cl0]
    hfrac.gene0 = names(which(frac0 > frac.thr))
    res = rbind(res, tf.diff.res0[feature %in% hfrac.gene0])
  }
  return(res)
} 

## plot enriched tf for chromvar, color by cluster
plot_enrich_tf0 <- function(sele.zscores, bc_clusters,
                           cluster.levels = NULL,
                           up.qt = 0.95, low.qt = 0.05,
                           ndownsample = 1000){
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample 
  set.seed(2019)
  bc_clusters.down = bc_clusters
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = as.character(bc_clusters$cluster), 
                          'barcode' = bc_clusters$barcode,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nc >= 3) color_cluster = getPalette(nc)
  if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
  names(color_cluster) = sort(unique(bc_clusters$cluster))
  
  ann_column = ann_column[order(ann_column$cluster), ]
  
  
  
  if(is.null(cluster.levels)) {
    
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = sort(unique(ann_column$cluster)))), ]
    
  }else{
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = cluster.levels)), ]
    
    color_cluster = color_cluster[cluster.levels]
  }
  ann_colors = list('cluster' = color_cluster)
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode = NULL
  
  ph <- pheatmap::pheatmap(sele.zscores, 
                           cluster_cols = F, cluster_rows = F, 
                           show_colnames = F, fontsize = 9,
                           annotation_col = ann_column, 
                           color = viridis(100),
                           annotation_colors = ann_colors, 
                           fontsize_row = 9)
  return(ph)
}

## virids style, downsample, match_cell, colomn color by sample/cluster
plot_enrich_tf1 <- function(sele.tfs, zscore.all, bc_clusters,
                           cluster.levels = NULL, 
                           cluster.col = NULL,
                           up.qt = 0.95, low.qt = 0.05,
                           ndownsample = 2000, match_cell = F,
                           reScale = F, cluster.rows = F,
                           color_style = 'virid'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$cluster))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(cl0 in unique(bc_clusters$cluster)){
      tmp = bc_clusters[cluster == cl0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(is.null(cluster.col)){
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
  }else{
    color_cluster = cluster.col
  }
  

  ## order 
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]
  if(is.null(cluster.levels)) {
    
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = sort(unique(ann_column$cluster)))), ]
    
  }else{
    ann_column = ann_column[order(factor(ann_column$cluster,
                                         levels = cluster.levels)), ]
    
    color_cluster = color_cluster[cluster.levels]
  }
  
  
  ann_colors = list('cluster' = color_cluster)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}



## virids style, downsample, match_cell, column color by sample & cluster
plot_enrich_tf2 <- function(sele.tfs, zscore.all, bc_clusters,
                           cluster.levels = NULL, sample.levels = NULL,
                           up.qt = 0.95, low.qt = 0.05,
                           ndownsample = 2000, 
                           match_cell = F, reScale = F, 
                           cluster.rows = F,
                           color_style = 'virid',
                           order.within = 'sample'){
  sele.zscores = zscore.all[sele.tfs, ]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$sample))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(sample0 in unique(bc_clusters$sample)){
      tmp = bc_clusters[sample == sample0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          'sample' = bc_clusters$sample,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  nc = length(unique(bc_clusters$cluster))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nc >= 3) color_cluster = getPalette(nc)
  if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
  names(color_cluster) = sort(unique(bc_clusters$cluster))

  nsample = length(unique(bc_clusters$sample))
  getPalette = colorRampPalette(brewer.pal(9, "Paired"))
  if(nsample < 3) color_sample = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nsample]
  if(nsample > 3) color_sample = getPalette(nsample)
  names(color_sample) = sample.levels
  
  ## order sample by age
  ann_column = ann_column[order(factor(ann_column$sample,
                                       levels = sample.levels[sample.levels %in% ann_column$sample])), ]
  ann_column = ann_column[order(factor(ann_column$cluster,
                                       levels = sort(unique(ann_column$cluster)))), ]

  if(order.within == 'sample'){
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
  }else{
    if(!is.null(sample.levels)){
      ann_column = ann_column[order(factor(ann_column$sample,
                                           levels = sample.levels)), ]
      
      color_sample = color_sample[sample.levels]
    }
    
    if(!is.null(cluster.levels)) {
      ann_column = ann_column[order(factor(ann_column$cluster,
                                           levels = cluster.levels)), ]
      
      color_cluster = color_cluster[cluster.levels]
    }
    
    
    
  }

  ann_colors = list('cluster' = color_cluster,
                    'sample' = color_sample)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
  
  
}

## virids style, downsample, match_cell, colomn color by a continous variable (time)
## and a cluster/sample
plot_enrich_tf3 <- function(sele.tfs, zscore.all, bc_clusters,
                            cluster.levels = NULL, cluster.col = NULL,
                            up.qt = 0.95, low.qt = 0.05,
                            ndownsample = 2000, 
                            match_cell = F, reScale = F, 
                            cluster.rows = F,
                            color_style = 'virid', 
                            order.withinCls = F){
  sele.zscores = zscore.all[sele.tfs, ]
  sele.zscores = sele.zscores[, bc_clusters$barcode]
  
  if(reScale) sele.zscores = t(scale(t(sele.zscores), center = T, scale = T))
  
  bc_clusters = data.table(bc_clusters)
  
  #downsample and match_cell
  ncell.cl = min(table(bc_clusters$cluster))
  set.seed(2020)
  bc_clusters.down = bc_clusters
  if(match_cell){
    bc_clusters.down = NULL
    for(cl0 in unique(bc_clusters$cluster)){
      tmp = bc_clusters[cluster == cl0]
      if(nrow(tmp) > ncell.cl) tmp = tmp[sort(sample((1:nrow(tmp)), ncell.cl)), ]
      bc_clusters.down = rbind(bc_clusters.down, tmp)
    }
  }
  
  if(!is.null(ndownsample) & ndownsample < nrow(bc_clusters.down)) 
    bc_clusters.down = bc_clusters.down[sort(sample((1:nrow(bc_clusters.down)), ndownsample)), ]
  
  bc_clusters = bc_clusters.down
  rr = bc_clusters$barcode[bc_clusters$barcode %in% colnames(sele.zscores)]
  sele.zscores = sele.zscores[, rr]
  
  ann_column = data.frame('cluster' = bc_clusters$cluster,
                          'barcode' = bc_clusters$barcode,
                          'pseudotime' = bc_clusters$pseudotime,
                          stringsAsFactors = F)
  rownames(ann_column) = bc_clusters$barcode
  
  up_cut = quantile(sele.zscores, up.qt, na.rm = T)
  low_cut = quantile(sele.zscores, low.qt, na.rm = T)
  sele.zscores[is.na(sele.zscores)] = 0
  low_cut = min(0, low_cut)
  sele.zscores[sele.zscores > up_cut] = up_cut
  sele.zscores[sele.zscores < low_cut] = low_cut
  
  if(is.null(cluster.col)){
    nc = length(unique(bc_clusters$cluster))
    getPalette = colorRampPalette(brewer.pal(9, "Paired"))
    if(nc >= 3) color_cluster = getPalette(nc)
    if(nc < 3) color_cluster = c("#A6CEE3", "#1F78B4", "#B2DF8A")[1:nc]
    names(color_cluster) = sort(unique(bc_clusters$cluster))
  }else{
    color_cluster = cluster.col
  }
  
  ann_column = ann_column[order(ann_column$pseudotime), ]
  
  if(order.withinCls & !is.null(cluster.levels)){
      color_cluster = color_cluster[cluster.levels]
      tmp = NULL
      for(cl0 in cluster.levels){
        ann_column0 = ann_column[ann_column$cluster == cl0, ]
        ann_column0 = ann_column0[order(ann_column0$pseudotime), ]
        tmp = rbind(tmp, ann_column0)
      }
      ann_column = tmp
      rownames(ann_column) = ann_column$barcode
      ann_column$cluster = factor(ann_column$cluster, levels = cluster.levels)
      
  }
  ## order by psudotime
  
  ann_colors = list('cluster' = color_cluster)
  
  sele.zscores = sele.zscores[, ann_column$barcode]
  ann_column$barcode <- NULL
  
  color_fun = viridis(100)
  if(color_style == 'purple-yellow') color_fun = PurpleAndYellow()
  
  ph <- pheatmap::pheatmap(sele.zscores, cluster_cols = F, 
                           cluster_rows = cluster.rows, 
                           show_colnames = F, fontsize = 10,
                           annotation_col = ann_column, 
                           color = color_fun,
                           annotation_colors = ann_colors, 
                           fontsize_row = 12)
  return(ph)
  
  
}


## prepare gene-by-group freq of expressed cells ####
seurat.rna.hd <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')
mtx = seurat.rna.hd@assays$RNA@counts
Ctypes = seurat.rna.hd$Ctype

efreq.ctype.hd <- sapply(unique(Ctypes), function(x) {
  
  cl_data <- mtx[, Ctypes == x]
  
  Matrix::rowMeans(cl_data > 0)
  
})

efreq.earlyprog.hd <- rowMeans(mtx[, Ctypes == 'HSPC' | Ctypes == 'LMPP'] > 0)

rm(seurat.rna.hd, mtx)

seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
mtx.mllr = seurat.rna.mllr[['RNA']]@counts
Ctype0 = seurat.rna.mllr$Ctype0
Ctype0[Ctype0 %in% c('NKT', 'T')] = 'T/NK'

pCtypes = seurat.rna.mllr$projCtype
efreq.ctype.mllr <- sapply(unique(pCtypes), function(x) {
  
  cl_data <- mtx.mllr[, Ctype0 == 'Blasts' & pCtypes == x]
  
  Matrix::rowMeans(cl_data > 0)
  
})

efreq.ctype0.mllr <- sapply(unique(Ctype0), function(x) {
  
  cl_data <- mtx.mllr[, Ctype0 == x]
  
  Matrix::rowMeans(cl_data > 0)
  
})


efreq.hspc1 = rowMeans(mtx.mllr[, Ctype0 == 'Progenitors'] > 0)
efreq.earlyprog.mllr = rowMeans(mtx.mllr[, Ctype0 == 'Blasts' & (pCtypes == 'HSPC-like' | pCtypes == 'LMPP-like')] > 0)

rm(seurat.rna.mllr, mtx.mllr)


## save data for future usage
save(efreq.ctype.hd, efreq.earlyprog.hd, file = 'MetaData/scRNA/efreq_ctypte_hd.RData')
save(efreq.ctype.mllr, efreq.hspc1, efreq.earlyprog.mllr,
     efreq.ctype0.mllr, file = 'MetaData/scRNA/efreq_ctypte_mllr.RData')


## prepare dscore and zscore ####
if(F){
  chrvar.all <- readRDS('chromVAR_Objects/chromVAR_MLLrwithPubAtac_usePeaks4TFcomp_allPeaks.rds')
  
  zscore.all = chrvar.all@assays@data$z
  dscore.all = chrvar.all@assays@data$deviations
   
  ## change to readable name
  rnames = rownames(zscore.all)
  nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
  nnames1 = sapply(rnames, function(x) unlist(strsplit(x, '_'))[1])
  rownames(zscore.all) = ifelse(grepl(nnames, pattern = 'LINE'),
                               nnames1, nnames)
  rownames(dscore.all) = rownames(zscore.all)
  
  ## ENSG00000187728 -- TCF24
  
  save(zscore.all, dscore.all, 
       file = 'chromVAR_Objects/zscore_dscore_all.RData')
  
  save(zscore.all, dscore.all, 
       file = 'chromVAR_Objects/zscore_dscore_all_allPeaks.RData')
}


## add metadata for each cell ####
load('chromVAR_Objects/zscore_dscore_all_allPeaks.RData')

seurat.hd <- readRDS('Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds')
seurat.mllr <- readRDS('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap10000.rds')
ann.hd = data.table(seurat.hd@meta.data, keep.rownames = T)
ann.mllr = data.table(seurat.mllr@meta.data, keep.rownames = T)
ann.hd = subset(ann.hd, select = c('rn', 'sample', 'seurat_ctype', 'pseudotime'))
ann.mllr = subset(ann.mllr, select = c('rn', 'sample', 'projCtype', 'projCtype_merge', 'Ctype0', 'pseudotime'))
ann.hd[, 'bc_name' := gsub('.1', '', rn, fixed = T), by = rn]
ann.hd[, 'bc_name' := gsub('.2', '', bc_name, fixed = T), by = rn]
ann.hd[, 'bc_name' := gsub('Pub3_', '', bc_name, fixed = T), by = rn]
ann.hd[, 'sample_cell' := paste0(sample, bc_name)]
ann.hd = subset(ann.hd, select = c('bc_name', 'sample', 'seurat_ctype', 'sample_cell',
                                   'pseudotime'))
names(ann.hd)[3] = 'Ctype'
fwrite(ann.hd, file = paste0('MetaData/scATAC/hd_bc_sample_inf.txt'),
      sep = '\t',  quote = F)

ann.mllr$bc_name = sapply(ann.mllr$rn, function(x) substr(x, nchar(x)-15, nchar(x)))
ann.mllr = subset(ann.mllr, select = c('bc_name', 'sample', 'projCtype', 'projCtype_merge', 'Ctype0',
                                       'pseudotime'))
ann.mllr[, 'sample_cell' := paste0(sample, bc_name)]
fwrite(ann.mllr, file = paste0('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt'),
       sep = '\t',  quote = F)


## compare blast projected cell type with corresponding normal cell types ####
load('chromVAR_Objects/zscore_dscore_all.RData')
ann.hd = fread('MetaData/scATAC/hd_bc_sample_inf.txt')
ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
ann.hspc1 = ann.mllr[Ctype0 == 'Progenitors']


b.types = c('HSPC', 'LMPP', 'CLP', 'Pre-pro-B', 'Pro-B', 'Pre-B', 'Immature-B', 
            'Mature-B')
b.types.like = paste0(b.types, '-like')
pairs.b.ctype = paste0(b.types, '_vs_', b.types.like)
pairs.other = c('HSPC1_vs_HSPC', 'EarlyProg_vs_EarlyProg-like')

res.diff = list()
for(pair0 in c(pairs.b.ctype, pairs.other)){
  if(pair0 %in% pairs.b.ctype){
    ctype0 = unlist(strsplit(pair0, split = '_vs_'))[1]
    ctype1 = paste0(ctype0, '-like')
    cls = rep(c(ctype0, ctype1),
              c(nrow(ann.hd[Ctype == ctype0]), nrow(ann.mllr[projCtype == ctype1 & Ctype0 == 'Blasts'])))
    cls = data.table('barcode' = c(ann.hd[Ctype == ctype0]$sample_cell, 
                                   ann.mllr[projCtype == ctype1 & Ctype0 == 'Blasts']$sample_cell),
                     'cluster' = cls)
  }
  if(pair0 == 'HSPC1_vs_HSPC'){
   
    cls = rep(c('HSPC', 'HSPC1'),
              c(nrow(ann.hd[Ctype == 'HSPC']), nrow(ann.hspc1)))
    cls = data.table('barcode' = c(ann.hd[Ctype == 'HSPC']$sample_cell, 
                                   ann.hspc1$sample_cell),
                     'cluster' = cls)
  }
  if(pair0 == 'HSPC1_vs_HSPC-like'){
    
    cls = rep(c('HSPC-like', 'HSPC1'),
              c(nrow(ann.mllr[Ctype0 == 'Blasts' & projCtype == 'HSPC-like']), 
                nrow(ann.hspc1)))
    cls = data.table('barcode' = c(ann.mllr[Ctype0 == 'Blasts' & projCtype == 'HSPC-like']$sample_cell, 
                                   ann.hspc1$sample_cell),
                     'cluster' = cls)
  }
  if(pair0 == 'HSPC1_vs_EarlyProg-like'){
    
    cls = rep(c('HSPC1', 'EarlyProg-like'),
              c(nrow(ann.hspc1), 
                nrow(ann.mllr[Ctype0 == 'Blasts' & projCtype %in% paste0(c('HSPC','LMPP'), '-like')])))
    cls = data.table('barcode' = c(ann.hspc1$sample_cell, 
                                   ann.mllr[Ctype0 == 'Blasts' & projCtype %in% 
                                              paste0(c('HSPC','LMPP'), '-like')]$sample_cell),
                     'cluster' = cls)
  }
  if(pair0 == 'HSPC1_vs_EarlyProg'){
    
    cls = rep(c('HSPC1', 'EarlyProg'),
              c(nrow(ann.hspc1), nrow(ann.hd[Ctype %in% c('HSPC', 'LMPP')])))
                
    cls = data.table('barcode' = c(ann.hspc1$sample_cell, 
                                   ann.hd[Ctype %in% c('HSPC', 'LMPP')]$sample_cell),
                     'cluster' = cls)
  }
  if(pair0 == 'EarlyProg_vs_EarlyProg-like'){
    
    cls = rep(c('EarlyProg', 'EarlyProg-like'),
              c(nrow(ann.hd[Ctype %in% c('HSPC', 'LMPP')]), 
                nrow(ann.mllr[Ctype0 == 'Blasts' & projCtype %in% paste0(c('HSPC','LMPP'), '-like')])))
    cls = data.table('barcode' = c(ann.hd[Ctype %in% c('HSPC', 'LMPP')]$sample_cell, 
                                   ann.mllr[Ctype0 == 'Blasts' & projCtype %in% 
                                            paste0(c('HSPC','LMPP'), '-like')]$sample_cell),
                     'cluster' = cls)
  }
  res.diff[[pair0]] <- runDiffMotifEnrich(dscore.all, clusters = cls, fdr = 0.05,
                                  max_cell_per_clust = 2000,
                                  topn = 200, min_frac_per_cluster = 0.2)
  rm(cls)
}

## filtered enriched list by expression 
load('MetaData/scRNA/efreq_ctypte_hd.RData')
load('MetaData/scRNA/efreq_ctypte_mllr.RData')
res.diff.filtered = list()
expr_cutoff = 0.2
for(pair0 in c(pairs.b.ctype, pairs.other)){
  res.diff0 = res.diff[[pair0]]
  sele.tfs0 = sapply(unique(res.diff0$cluster1), function(x){
    res.diff0 = res.diff0[order(-mean1)]
    tfs = res.diff0[cluster1 == x]$feature
    if(x %in% colnames(efreq.ctype.mllr)) efreq0 = efreq.ctype.mllr[, x]
    if(x %in% colnames(efreq.ctype.hd)) efreq0 = efreq.ctype.hd[, x]
    if(x == 'HSPC1') efreq0 = efreq.hspc1
    if(x == 'EarlyProg') efreq0 = efreq.earlyprog.hd
    if(x == 'EarlyProg-like') efreq0 = efreq.earlyprog.mllr
    ## filter by expr frac
    tfs = tfs[tfs %in% names(efreq0)]
    tfs = tfs[efreq0[tfs] > expr_cutoff]
    
    return(tfs)
  })
  if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('c', sele.tfs0)
  res.diff.filtered[[pair0]] = res.diff0[feature %in% sele.tfs0]
}

if(F){
  ## use less stringent rule to check more on the list
  res.diff.filtered = list()
  expr_cutoff = 0.05
  for(pair0 in c(pairs.b.ctype, pairs.other)){
    res.diff0 = res.diff[[pair0]]
    sele.tfs0 = sapply(unique(res.diff0$cluster1), function(x){
      res.diff0 = res.diff0[order(-mean1)]
      tfs = res.diff0[cluster1 == x]$feature
      if(x %in% colnames(efreq.ctype.mllr)) efreq0 = efreq.ctype.mllr[, x]
      if(x %in% colnames(efreq.ctype.hd)) efreq0 = efreq.ctype.hd[, x]
      if(x == 'HSPC1') efreq0 = efreq.hspc1
      if(x == 'EarlyProg') efreq0 = efreq.earlyprog.hd
      if(x == 'EarlyProg-like') efreq0 = efreq.earlyprog.mllr
      ## filter by expr frac
      tfs = tfs[tfs %in% names(efreq0)]
      tfs = tfs[efreq0[tfs] > expr_cutoff]
      
      return(tfs)
    })
    if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('c', sele.tfs0)
    res.diff.filtered[[pair0]] = res.diff0[feature %in% sele.tfs0]
  }
  #output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
  #                    expr_cutoff)
  dir.create(output_dir, showWarnings = F)
  blasts.diff = res.diff.filtered[1:8]
  blasts.diff.up = lapply(blasts.diff, function(x) x[grepl(cluster1, pattern = 'like')])
  blasts.diff.down = lapply(blasts.diff, function(x) x[!grepl(cluster1, pattern = 'like')])
  blasts.diff.up = do.call('rbind', blasts.diff.up)
  blasts.diff.down = do.call('rbind', blasts.diff.down)
  
}
## commonly enriched in blast cells
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)
dir.create(output_dir, showWarnings = F)
blasts.diff = res.diff.filtered[1:8]
blasts.diff.up = lapply(blasts.diff, function(x) x[grepl(cluster1, pattern = 'like')])
blasts.diff.down = lapply(blasts.diff, function(x) x[!grepl(cluster1, pattern = 'like')])
blasts.diff.up = do.call('rbind', blasts.diff.up)
blasts.diff.down = do.call('rbind', blasts.diff.down)
up.tfs = sort(table(blasts.diff.up$feature), decreasing = T)
down.tfs = sort(table(blasts.diff.down$feature), decreasing = T)

# plot a barplot for top TFs
if(T){
  tf_freq.up = data.table('TF' = names(up.tfs),
                       'N' = as.numeric(up.tfs))
  tf_freq.down = data.table('TF' = names(down.tfs),
                       'N' = as.numeric(down.tfs))
  
  tf_freq.up = tf_freq.up[order(-N), ]
  tf_freq.up = tf_freq.up[N >= 3]
  tf_freq.up = tf_freq.up[order(N), ]
  tf_freq.up$TF = factor(tf_freq.up$TF, levels = tf_freq.up$TF)
  
  p_up <- ggplot(tf_freq.up, aes(y = N, x = TF)) +
    geom_bar(width = 0.7, stat = 'identity', fill = '#1F78B4') +
    ggtitle("Commonly Enriched TF in Blasts") + theme_classic() + 
    theme(legend.position = 'bottom', legend.direction = "horizontal") + 
    coord_flip()  +
    xlab('') + ylab('# of Stages')
  

  tf_freq.down = tf_freq.down[N >= 3]
  tf_freq.down = tf_freq.down[order(N), ]
  tf_freq.down$TF = factor(tf_freq.down$TF, levels = tf_freq.down$TF)
  
  p_down <- ggplot(tf_freq.down, aes(y = N, x = TF)) +
    geom_bar(width = 0.7, stat = 'identity', fill = '#1F78B4') +
    ggtitle("Commonly Depleted TF in Blasts") + theme_classic() + 
    theme(legend.position = 'bottom', legend.direction = "horizontal") + 
    coord_flip()  +
    xlab('') + ylab('# of Stages')
  
  ggsave(p_up, filename = paste0(output_dir, '/Commonly_Enriched_TFs_InBlasts.eps'),
         width = 6, height = 7, device = 'eps')
  ggsave(p_down, filename = paste0(output_dir, '/Commonly_Depleted_TFs_InBlasts.eps'),
         width = 6, height = 10, device = 'eps')
  saveRDS(tf_freq.up, file = paste0(output_dir, '/commonly_enriched_tfs.rds'))
  saveRDS(tf_freq.down, file = paste0(output_dir, '/commonly_depleted_tfs.rds'))
}


## plot TFs heatmap commonly enriched in blast vs normal ####
up.tfs = names(which(up.tfs > 3))
down.tfs = names(which(down.tfs > 3))

sele.tfs = c(up.tfs, down.tfs)

cls = data.table('barcode' = c(ann.hd$sample_cell, ann.mllr$sample_cell),
                 'cluster' = c(ann.hd$Ctype, ann.mllr$projCtype),
                 'sample' = rep(c('HD', 'Blasts'), c(nrow(ann.hd), nrow(ann.mllr))))
ordered.clusters = c(b.types[1:8], paste0(b.types[1:8], '-like'))
ordered.sample = c('HD', 'Blasts')
cls = cls[cluster %in% ordered.clusters ]
ph <- plot_enrich_tf2(sele.tfs, zscore.all, bc_clusters = cls,
                     up.qt = 0.96,
                     ndownsample = 1000,
                     cluster.levels = ordered.clusters, cluster.rows = T,
                     sample.levels = ordered.sample, reScale = T,
                     match_cell = T, color_style = 'purple-yellow')

#ggsave(ph, filename = paste0(output_dir, '/TF_enrichedIn4Stages_blastVSnormal_purpleyellow.eps'),
#       width = 12, height = 12, device = 'eps')
ggsave(ph, filename = paste0(output_dir, '/TF_enrichedIn4Stages_blastVSnormal_purpleyellow_v2.eps'),
       width = 12, height = 12, device = 'eps')

## plot HSPC vs HSPC1 ####
hspc.diff = res.diff.filtered[['HSPC1_vs_HSPC']]
setkey(hspc.diff, cluster1, pv_adjust)
cls = data.table('barcode' = c(ann.hd[Ctype == 'HSPC']$sample_cell, 
                               ann.hspc1$sample_cell),
                 'cluster' = c(ann.hd[Ctype == 'HSPC']$Ctype, 
                               rep('HSPC1', nrow(ann.hspc1))),
                  'sample' = c(rep('HD', nrow(ann.hd[Ctype == 'HSPC'])),
                               ann.hspc1$sample))
ordered.sample = c('HD', 'MLLr876533', 'MLLr882304', 'MLLr875706',
                   'MLLr1154', 'MLLr870684', 'MLLr879583',
                   'MLLr876545', 'MLLr874013', 'MLLr879440',
                   'MLLr878289', 'MLLr871427', 'MLLr875703',
                   'MLLr877476', 'MLLr881823', 'MLLr877780',
                   'MLLr879339', 'MLLr878501', 'MLLr878516')
ordered.clusters = c('HSPC', 'HSPC1')
hspc.diff[, 'delta' := -1*(mean1 - mean0)/mean1]
setkey(hspc.diff, cluster1, delta)
hspc.diff$delta = -hspc.diff$delta
hspc.diff = hspc.diff[mean1 > 0.012]

ph.py <- plot_enrich_tf1(hspc.diff$feature, zscore.all, bc_clusters = cls,
                      up.qt = 0.97,
                      ndownsample = 1000, match_cell = T,
                      cluster.levels = ordered.clusters, reScale = T,
                      color_style = 'purple-yellow')
ggsave(ph.py, filename = paste0(output_dir, '/enriched_TFs_HSPCvsHSPC1_purpleyellow.eps'),
                width = 12, height = 12, device = 'eps')


## comparing blast with normal cells in patients ####
## prepare dscore and zscore
if(F){
  
  chrvar.mllr <- readRDS('chromVAR_Objects/chromVAR_pool_18MLLr_TFIDF_allPeaks.rds')
  zscore.mllr = chrvar.mllr@assays@data$z
  dscore.mllr = chrvar.mllr@assays@data$deviations
  
  ## change to readable name
  rnames = rownames(zscore.mllr)
  nnames = sapply(rnames, function(x) unlist(strsplit(x, '_'))[3])
  nnames1 = sapply(rnames, function(x) unlist(strsplit(x, '_'))[1])
  rownames(zscore.mllr) = ifelse(grepl(nnames, pattern = 'LINE'),
                                  nnames1, nnames)
  rownames(dscore.mllr) = rownames(zscore.mllr)
  
  ## ENSG00000187728 -- TCF24
  
  #save(zscore.mllr, dscore.mllr, 
  #     file = 'chromVAR_Objects/zscore_dscore_18MLLr.RData')
  save(zscore.mllr, dscore.mllr, 
       file = 'chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')
  
}
load('chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')

ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
load('MetaData/scRNA/efreq_ctypte_mllr.RData')
ctype0.count = table(ann.mllr$Ctype0)

cls = subset(ann.mllr, select = c(Ctype0, sample_cell, sample))
names(cls) = c('cluster', 'barcode', 'sample')
dscore.mllr[is.na(dscore.mllr)] = 0

## change cell names to be consistent
cnames = colnames(dscore.mllr)
cnames = sapply(cnames, function(x) gsub('_scATAC', '', x))
cnames = sapply(cnames, function(x) gsub('MLL_', 'MLLr', x))
all(cnames %in% cls$barcode)
colnames(dscore.mllr) = cnames
colnames(zscore.mllr) = cnames
cls = cls[barcode %in% cnames]
tf.mllr <- runDiffMotifEnrich(dscore.mllr, clusters = cls, fdr = 0.05,
                               max_cell_per_clust = 2000,
                               topn = 200, min_frac_per_cluster = 0.2)

if(F){
  ## relax TFs for GRN
  expr_cutoff = 0.1
  topn = 20
  tf.mllr.sele = tf.mllr[cluster1 %in% c('Blasts', 'Progenitors', 'Mature_B',
                                         'T/NK', 'Monocytes')]
  sele.tfs0 = lapply(c('Blasts', 'Progenitors', 'Mature_B',
                       'T/NK', 'Monocytes'), function(x){
    res.diff0 = tf.mllr.sele[order(-mean1)]
    res.diff0 = res.diff0[cluster1 == x]
    tfs = res.diff0[cluster1 == x]$feature
    efreq0 = efreq.ctype0.mllr[, x]
    ## filter by expr frac
    tfs = tfs[tfs %in% names(efreq0)]
    tfs = tfs[efreq0[tfs] > expr_cutoff]
    tfs = res.diff0[1:topn, ]
    return(tfs)
  })
  names(sele.tfs0) = c('Blasts', 'Progenitors', 'Mature_B',
                       'T/NK', 'Monocytes')
  saveRDS(sele.tfs0, file = paste0('EP_Prediction/Ctype0_enriched_tfs_exprfrac0.1.rds'))
  
}
expr_cutoff = 0.2
tf.mllr.sele = tf.mllr[cluster1 %in% c('Blasts')]
sele.tfs0 = sapply(c('Blasts'), function(x){
  res.diff0 = tf.mllr.sele[order(-mean1)]
  tfs = res.diff0[cluster1 == x]$feature
  efreq0 = efreq.ctype0.mllr[, x]
  ## filter by expr frac
  tfs = tfs[tfs %in% names(efreq0)]
  tfs = tfs[efreq0[tfs] > expr_cutoff]
  
  return(tfs)
})

if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('c', sele.tfs0)
sele.tfs0 = unique(sele.tfs0)

tf.mllr.sele = tf.mllr.sele[feature %in% sele.tfs0]
tf.mllr.sele[, 'delta' := (mean1 - mean0)/mean1]
tf.mllr.sele = tf.mllr.sele[order(-delta), ]
cls = cls[cluster != 'pDC']


ordered.sample = c('MLLr876533', 'MLLr882304', 'MLLr875706',
                   'MLLr1154', 'MLLr870684', 'MLLr879583',
                   'MLLr876545', 'MLLr874013', 'MLLr879440',
                   'MLLr878289', 'MLLr871427', 'MLLr875703',
                   'MLLr877476', 'MLLr881823', 'MLLr877780',
                   'MLLr879339', 'MLLr878501', 'MLLr878516')

ph <- plot_enrich_tf2(tf.mllr.sele$feature, zscore.mllr,
                      bc_clusters = cls,
                      up.qt = 0.96,
                      ndownsample = 1000, match_cell = F,
                      cluster.levels = c('Progenitors', 'Blasts', 
                                         'Mature_B', 'Monocytes',
                                         'T/NK'),
                      sample.levels = ordered.sample,
                      reScale = T, cluster.rows = T, 
                      order.within = 'cluster',
                      color_style = 'purple-yellow')
ggsave(ph, filename = paste0(output_dir, '/enriched_TFs_Blast_Normal.eps'),
       width = 12, height = 12, device = 'eps')

## show top TFs for all ctype0 ####

ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')


load('MetaData/scRNA/efreq_ctypte_mllr.RData')
ctype0.count = table(ann.mllr$Ctype0)

expr_cutoff = 0.2
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)
tf.mllr.sele = tf.mllr
topn = 15

## plot evenly across cell types 
sele.tfs0 = lapply(c('Progenitors', 'Blasts', 
                     'Mature_B', 'Monocytes',
                     'T/NK'), function(x){
  res.diff0 = tf.mllr.sele[cluster1 == x]
  res.diff0 = res.diff0[order(-mean1)]
  efreq0 = efreq.ctype0.mllr[, x]
  ## filter by expr frac
  res.diff0 = res.diff0[feature %in% names(efreq0)]
  tfs = res.diff0$feature
  tfs.expr = names(which(efreq0[tfs] > expr_cutoff))
  res.diff0 = res.diff0[feature %in% tfs.expr]
  if(nrow(res.diff0) > topn) res.diff0 = res.diff0[1:topn, ]
  return(res.diff0)
})
names(sele.tfs0) = c('Progenitors', 'Blasts', 
                     'Mature_B', 'Monocytes',
                     'T/NK')
saveRDS(sele.tfs0, file = paste0(output_dir, '/Ctype0_top_enriched_tfs.rds'))

if(class(sele.tfs0) == 'list') sele.tfs = do.call('rbind', sele.tfs0)
sele.tfs = sele.tfs$feature

sele.tfs = unique(c(sele.tfs0[['Progenitors']]$feature, sele.tfs0[['Blasts']]$feature,
                    sele.tfs0[['Mature_B']]$feature, sele.tfs0[['Monocytes']]$feature))
sele.tfs = c(setdiff(sele.tfs, c('ELF1', 'ETS1', 'ELK4')), sele.tfs0[['T/NK']]$feature)

cls_color <- c(  "Blasts" = "#E3191C", "Monocytes" = "#347EB5", 
  "Mature_B" = "#4CAC47",  "T/NK" = "#A35728", "Progenitors" = "#924C9F")

ph <- plot_enrich_tf1(unique(sele.tfs), zscore.mllr,
                      bc_clusters = cls,
                      up.qt = 0.96, low.qt = 0.1,
                      ndownsample = 1000, match_cell = T,
                      cluster.levels = c('Progenitors', 'Blasts', 
                                         'Mature_B', 'Monocytes',
                                         'T/NK'),
                      cluster.col = cls_color,
                      reScale = T, cluster.rows = F, 
                      color_style = 'purple-yellow')
ggsave(ph, filename = paste0(output_dir, '/enriched_TFs_Ctype0.eps'),
       width = 14, height = 14, device = 'eps')

## plot use real cell type composition
sele.tfs = unique(c(sele.tfs0[['Blasts']]$feature, sele.tfs0[['T/NK']]$feature,
                 sele.tfs0[['Mature_B']]$feature, sele.tfs0[['Monocytes']]$feature,
                 sele.tfs0[['Progenitors']]$feature))


ph1 <- plot_enrich_tf1(unique(sele.tfs), zscore.mllr,
                      bc_clusters = cls,
                      up.qt = 0.96,
                      ndownsample = 1000, match_cell = F,
                      cluster.levels = c( 'Blasts',  'T/NK',
                                         'Mature_B', 'Monocytes',
                                         'Progenitors'),
                      cluster.col = cls_color,
                      reScale = T, cluster.rows = F, 
                      color_style = 'purple-yellow')
ggsave(ph1, filename = paste0(output_dir, '/enriched_TFs_Ctype0_unevenComposition.eps'),
       width = 14, height = 14, device = 'eps')




## comparing blasts within patients among merged projected cell type ####

## prepare Early-prog -like expression fraction
## Early-prog-like: HSPC-like, LMPP-like, CLP-like
seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
mtx.mllr = seurat.rna.mllr[['RNA']]@counts
tmp = subset(seurat.rna.mllr, Ctype_Stage %in% c('HSPC-like', 'LMPP-like', 'CLP-like'))
efreq.earlyprog.mllr1 = rowMeans(tmp@assays$RNA@counts > 0)

load('chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')

ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
load('MetaData/scRNA/efreq_ctypte_mllr.RData')
ann.mllr = ann.mllr[Ctype0 == 'Blasts']
b.types = paste0(c('HSPC', 'LMPP', 'CLP', 'Pre-pro-B', 'Pre-B', 'Pro-B', 
            'Immature-B', 'Mature-B'), '-like')
ann.mllr = ann.mllr[projCtype %in% b.types]

ann.mllr[, 'mergedCtype' := ifelse(projCtype %in% c('HSPC-like', 'LMPP-like', 'CLP-like'),
                                   'Early-Prog-like', projCtype)]

cls = subset(ann.mllr, select = c(mergedCtype, sample_cell, pseudotime))
names(cls) = c('cluster', 'barcode', 'pseudotime')
cls = cls[order(pseudotime)]
dscore.mllr[is.na(dscore.mllr)] = 0

## change cell names to be consistent
cnames = colnames(dscore.mllr)
cnames = sapply(cnames, function(x) gsub('_scATAC', '', x))
cnames = sapply(cnames, function(x) gsub('MLL_', 'MLLr', x))
colnames(dscore.mllr) = cnames
colnames(zscore.mllr) = cnames
cls = cls[barcode %in% cnames]
tf.mllr <- runDiffMotifEnrich(dscore.mllr, clusters = cls, fdr = 0.05,
                              max_cell_per_clust = 2000,
                              topn = 200, min_frac_per_cluster = 0.2)

expr_cutoff = 0.2
topn = 15
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)
sele.tfs0 = lapply(c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                     'Pre-B-like', 'Immature-B-like', 'Mature-B-like'), function(x){
                       res.diff0 = tf.mllr[cluster1 == x]
                       res.diff0 = res.diff0[order(-mean1)]
                       if(x == 'Early-Prog-like'){
                         efreq0 = efreq.earlyprog.mllr1
                       }else{
                         efreq0 = efreq.ctype.mllr[, x]
                       }
                       ## filter by expr frac
                       res.diff0 = res.diff0[feature %in% names(efreq0)]
                       tfs = res.diff0$feature
                       tfs.expr = names(which(efreq0[tfs] > expr_cutoff))
                       res.diff0 = res.diff0[feature %in% tfs.expr]
                       if(nrow(res.diff0) > topn) res.diff0 = res.diff0[1:topn, ]
                       return(res.diff0)
                     })
names(sele.tfs0) = c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                     'Pre-B-like', 'Immature-B-like', 'Mature-B-like')


if(F){
  ## select less stringent TF for constructing GRN
  expr_cutoff = 0.1
  topn = 20
  output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                      expr_cutoff)
  sele.tfs0 = lapply(c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                       'Pre-B-like', 'Immature-B-like', 'Mature-B-like'), function(x){
                         res.diff0 = tf.mllr[cluster1 == x]
                         res.diff0 = res.diff0[order(-mean1)]
                         if(x == 'Early-Prog-like'){
                           efreq0 = efreq.earlyprog.mllr1
                         }else{
                           efreq0 = efreq.ctype.mllr[, x]
                         }
                         ## filter by expr frac
                         res.diff0 = res.diff0[feature %in% names(efreq0)]
                         tfs = res.diff0$feature
                         tfs.expr = names(which(efreq0[tfs] > expr_cutoff))
                         res.diff0 = res.diff0[feature %in% tfs.expr]
                         if(nrow(res.diff0) > topn) res.diff0 = res.diff0[1:topn, ]
                         return(res.diff0)
                       })
  names(sele.tfs0) = c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                       'Pre-B-like', 'Immature-B-like', 'Mature-B-like')
  
  saveRDS(sele.tfs0, file = paste0('EP_Prediction/Blasts_MergeredStage_enriched_tfs_exprfrac0.1.rds'))
  
  
}

if(class(sele.tfs0) == 'list') sele.tfs = do.call('rbind', sele.tfs0)
sele.tfs = sele.tfs$feature
ctype_color <- c( 
  "Early-Prog-like" = "#E3191C", "Pre-pro-B-like" = "#347EB5",
  "Pro-B-like" = "#4CAC47", "Pre-B-like" = "#A35728", 
  "Immature-B-like" = "#924C9F", "Mature-B-like" = "979797"
)
ph0 <- plot_enrich_tf3(unique(sele.tfs), zscore.mllr,
                       bc_clusters = cls,
                       up.qt = 0.92, low.qt = 0.1,
                       ndownsample = 1000, match_cell = T,
                       cluster.levels = c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                                          'Pre-B-like', 'Immature-B-like', 'Mature-B-like'),
                       cluster.col = ctype_color,
                       reScale = T, cluster.rows = F, 
                       color_style = 'purple-yellow', order.withinCls = T)
ggsave(ph0, filename = paste0(output_dir, '/enriched_TFs_mergedStage_withinPatients.eps'),
       width = 14, height = 14, device = 'eps')


## regroup by projected lineage ####
tf.enrich.list = list()
fig.dir = paste0('Figures/scATAC/TF_Enrich/BySample_ByLineage')
dir.create(fig.dir, showWarnings = F)
b.types = c('HSPC', 'LMPP', 'CLP', 'Pre-pro-B', 'Pre-B', 'Pro-B', 
            'Immature-B', 'Mature-B')
dc.types = c('pDC', 'cDC', 'DC_Progenitor')
m.types = c('Mono', 'GMP')

for(sample0 in sampleNames){
  ## load dscore and zscore
  load(paste0('chromVAR_Objects/zscore_dscore_', sample0, '.RData'))
  cls = subset(mllr.ann, sample_patient==sample0, 
               select = c("patient_bc", "ctype_patient", "ctype_healthy0"))
  cls = cls[ctype_patient == 'Blast']
  cls$cluster = 'Other'
  cls$cluster = ifelse(cls$ctype_healthy0 %in% b.types, 'B-Lineage', cls$cluster)
  cls$cluster = ifelse(cls$ctype_healthy0 %in% m.types, 'M-Lineage', cls$cluster)
  cls$cluster = ifelse(cls$ctype_healthy0 %in% dc.types, 'DC-Lineage', cls$cluster)
  cls = subset(cls, select = c('patient_bc', 'cluster'))
  names(cls)[1] = 'barcode'
  tf.diff0 = runDiffMotifEnrich(dscore.mllr, clusters = cls,
                                topn = 200, 
                                min_frac_per_cluster = 0.1,
                                max_cell_per_clust = 1000)
  

  
  tf.diff0$sample = sample0
  tf.enrich.list[[sample0]] = tf.diff0
  
  ## plot enriched TFs in M-lineage per sample
  tf0 = tf.diff0[cluster1 == 'M-Lineage']
  if(nrow(tf0) == 0) next
  sele.zscore = zscore.mllr[tf0$feature, ]
  ph <- plot_enrich_tf(sele.zscore, bc_clusters = cls,
                       ndownsample = 1000)
  hh = floor(nrow(sele.zscore)/6)
  
  pfname = paste0(fig.dir, '/Heatmap_EnrichedTF_blasts_',
                  sample0, '_diffBlasts.eps')
  ggsave(ph, filename = pfname, device = 'eps', 
         height = max(4, hh),
         width = 9)
}

tf.enrich.sum = do.call('rbind', tf.enrich.list)
tf.enrich.sum = tf.enrich.sum[cluster1 == 'M-Lineage']
tf.freq <- sort(table(tf.enrich.sum$feature), decreasing = T)
tf.freq.sele = data.table(tf.freq[tf.freq >= 6])
names(tf.freq.sele)[1] = 'TF'
fwrite(tf.freq.sele, sep = '\t',
       file = paste0(fig.dir, '/topFreqTF_M_Lineage_blasts.csv'))


## pool all MLLr patients
dscore.comb <- readRDS("chromVAR_Objects/pool_12MLLr_chromVAR_dev.rds")
zscore.comb <- readRDS("chromVAR_Objects/pool_12MLLr_chromVAR_zscore.rds")
new.bcs =  sapply(colnames(dscore.comb), function(x){
  x <- gsub('CPTCA_Pilot', '1154', x)
  paste0('MLLr', x)
}) ## rename barcodes to avoid duplicates
names(new.bcs) <- NULL
colnames(zscore.comb) = colnames(dscore.comb) = new.bcs

mllr.ann[, 'new_patient_bc' := paste0(sample_patient, patient_bc)]
cls = subset(mllr.ann, select = c("new_patient_bc", "ctype_patient", "ctype_healthy0"))
cls = cls[ctype_patient == 'Blast']
cls$cluster = 'Other'
cls$cluster = ifelse(cls$ctype_healthy0 %in% b.types, 'B-Lineage', cls$cluster)
cls$cluster = ifelse(cls$ctype_healthy0 %in% m.types, 'M-Lineage', cls$cluster)
cls$cluster = ifelse(cls$ctype_healthy0 %in% dc.types, 'DC-Lineage', cls$cluster)
cls = subset(cls, select = c('new_patient_bc', 'cluster'))
names(cls)[1] = 'barcode'
tf.diff0 = runDiffMotifEnrich(dscore.comb, clusters = cls,
                              topn = 200, 
                              min_frac_per_cluster = 0.1,
                              max_cell_per_clust = 1000)
m.tf <- tf.diff0[cluster1 == 'M-Lineage']

ph <- plot_enrich_tf(zscore.comb[m.tf$feature[1:50], ], 
                     bc_clusters = cls, up.qt = 0.95,
               ndownsample = 1000)

ggsave(ph, filename = 'Figures/scATAC/TF_Enrich/topEnrichedTFs_M_Lineage_blasts.eps', 
       device = 'eps', 
       height = 7, width = 9)


## comparing young vs old blast in patients ####
## prepare dscore and zscore
load('chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')
ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna.mllr$projCtype[seurat.rna.mllr$projCtype %in% c('HSPC-like', 'LMPP-like',
                                                           'CLP-like')] <- 'Early-Prog-like'

seurat.rna.mllr$Age_Group = ifelse(seurat.rna.mllr$Age %in% c('10M', '11M', '7M',
                                                              '8M', '9M'), 'Old', 'Young')
seurat.rna.mllr = subset(seurat.rna.mllr, Ctype0 == 'Blasts')

## change cell names to be consistent
cnames = colnames(dscore.mllr)
cnames = sapply(cnames, function(x) gsub('_scATAC', '', x))
cnames = sapply(cnames, function(x) gsub('MLL_', 'MLLr', x))

colnames(dscore.mllr) = cnames
colnames(zscore.mllr) = cnames

younger.g <-  c('MLLr876533', 'MLLr882304', 'MLLr875706',
                'MLLr1154', 'MLLr870684', 'MLLr879583',
                'MLLr876545', 'MLLr874013', 'MLLr879440',
                'MLLr878289', 'MLLr871427')
ann.mllr$Age_Group = ifelse(ann.mllr$sample %in% younger.g, 'Young', 'Old')
ann.mllr$projCtype[ann.mllr$projCtype %in% c('HSPC-like', 'LMPP-like',
                                                           'CLP-like')] <- 'Early-Prog-like'

expr_cutoff = 0.1 
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)
ph = ph_bulk = list()
for(stage_type in c('pool', 'Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                    'Pre-B-like', 'Immature-B-like', 'Mature-B-like')){
  # prepare expression freq
  seurat.rna = seurat.rna.mllr
  if(stage_type != 'pool') seurat.rna = subset(seurat.rna, projCtype == stage_type)
  
  mtx.mllr = seurat.rna[['RNA']]@counts
  
  efreq.age.mllr <- sapply(c('Young', 'Old'), function(x) {
    
    cl_data <- mtx.mllr[, seurat.rna$Age_Group == x]
    
    Matrix::rowMeans(cl_data > 0)
    
  })
  if(stage_type == 'pool') save(efreq.age.mllr, file = 'MetaData/scRNA/efreq_ctype_YoungvsOld.RData')
  
  ann.mllr = ann.mllr[Ctype0 == 'Blasts']
  cls = subset(ann.mllr, select = c(Age_Group, sample_cell, sample, pseudotime))
  if(stage_type != 'pool')  cls = subset(ann.mllr, projCtype == stage_type,
                                         select = c(Age_Group, sample_cell, sample, pseudotime))
  names(cls) = c('cluster', 'barcode', 'sample', 'pseudotime')
  dscore.mllr[is.na(dscore.mllr)] = 0
  
 
  cls = cls[barcode %in% cnames]
  mcell = 2000
  if(stage_type %in% c('Early-Prog-like', 'Immature-B-like', 'Mature-B-like'))
    mcell = 100
  
  tf.mllr <- runDiffMotifEnrich(dscore.mllr, clusters = cls, fdr = 0.05,
                                max_cell_per_clust = mcell,
                                topn = 200, min_frac_per_cluster = 0.2)
  
  
  tf.mllr.sele = tf.mllr
  sele.tfs0 = lapply(c('Young', 'Old'), function(x){
    res.diff0 = tf.mllr.sele[order(-mean1)]
    tfs = res.diff0[cluster1 == x]$feature
    efreq0 = efreq.age.mllr[, x]
    ## filter by expr frac
    tfs = tfs[tfs %in% names(efreq0)]
    tfs = tfs[efreq0[tfs] > expr_cutoff]
    
    return(res.diff0[feature %in% tfs, ])
  })
  names(sele.tfs0) = c('Young', 'Old')
  if(stage_type == 'pool') saveRDS(sele.tfs0, "EP_Prediction/Blasts_YoungvsOld_enriched_tfs_exprfrac0.1.rds")
  
  if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('rbind', sele.tfs0)
  
  tf.mllr.sele = tf.mllr.sele[feature %in% sele.tfs0$feature]
  tf.mllr.sele = tf.mllr.sele[, delta := (mean1 - mean0)]
  tf.mllr.sele = tf.mllr.sele[order(-delta), ]
  nshow = min(50, nrow(tf.mllr.sele))
  tf.mllr.sele = tf.mllr.sele[1:nshow, ]
  
  age_color <- c( 
    'Young' = '#A6CEE3',
    'Old' = '#1F78B4'
  )
  ndownsample = min(1000, nrow(cls))
  ph[[stage_type]] <- plot_enrich_tf3(unique(tf.mllr.sele$feature), zscore.mllr,
                         bc_clusters = cls,
                         up.qt = 0.925, low.qt = 0.15,
                         ndownsample = ndownsample, match_cell = T,
                         cluster.levels = c('Young', 'Old'),
                         cluster.col = age_color,
                         reScale = T, cluster.rows = T, 
                         color_style = 'purple-yellow', 
                         order.withinCls = T)
  hmfile = paste0(output_dir, '/enriched_TFs_YoungVsOld.eps')
  if(stage_type != 'pool') hmfile = paste0(output_dir, '/enriched_TFs_YoungVsOld_', stage_type, '.eps')
  ht = ceiling(12 * nrow(tf.mllr.sele)/ 50 )
  ggsave(ph[[stage_type]], filename = hmfile,
         width = 12, height = 12, device = 'eps')
  
  ## plot group mean
  zscore.mllr1 = zscore.mllr[, cls$barcode]
  zscore.mllr1 = t(scale(t(zscore.mllr), center = T, scale = T))
    
  avg.y <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, cls[cluster=='Young']$barcode])
  avg.o <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, cls[cluster=='Old']$barcode])
  
  avg.score = cbind(avg.y, avg.o)
  colnames(avg.score) = c('Young', 'Old')
  ph_bulk[[stage_type]] <- pheatmap::pheatmap(avg.score, cluster_cols = F, 
                                              cluster_rows = T)
  
  hmfile_avg = paste0(output_dir, '/enriched_TFs_YoungVsOld_avg.eps')
  if(stage_type != 'pool') hmfile_avg = paste0(output_dir, '/enriched_TFs_YoungVsOld_', stage_type, '_avg.eps')
  ggsave(ph_bulk[[stage_type]], filename = hmfile_avg,
         width = 4, height = ht, device = 'eps')
  
}


## update young vs old blast (by stage and plot together in pseudo-bulk level) ####
## prepare dscore and zscore
load('chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')
ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna.mllr$projCtype[seurat.rna.mllr$projCtype %in% c('HSPC-like', 'LMPP-like',
                                                           'CLP-like')] <- 'Early-Prog-like'

seurat.rna.mllr$Age_Group = ifelse(seurat.rna.mllr$Age %in% c('10M', '11M', '7M',
                                                              '8M', '9M'), 'Old', 'Young')
seurat.rna.mllr = subset(seurat.rna.mllr, Ctype0 == 'Blasts')

## change cell names to be consistent
cnames = colnames(dscore.mllr)
cnames = sapply(cnames, function(x) gsub('_scATAC', '', x))
cnames = sapply(cnames, function(x) gsub('MLL_', 'MLLr', x))

colnames(dscore.mllr) = cnames
colnames(zscore.mllr) = cnames

younger.g <-  c('MLLr876533', 'MLLr882304', 'MLLr875706',
                'MLLr1154', 'MLLr870684', 'MLLr879583',
                'MLLr876545', 'MLLr874013', 'MLLr879440',
                'MLLr878289', 'MLLr871427')
ann.mllr$Age_Group = ifelse(ann.mllr$sample %in% younger.g, 'Young', 'Old')
ann.mllr$projCtype[ann.mllr$projCtype %in% c('HSPC-like', 'LMPP-like',
                                             'CLP-like')] <- 'Early-Prog-like'

expr_cutoff = 0.1
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)
dir.create(output_dir, showWarnings = F)
tf.res.list = avg.score.list = list()
for(stage_type in c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
                    'Pre-B-like', 'Immature-B-like', 'Mature-B-like')){
  # prepare expression freq
  seurat.rna = seurat.rna.mllr
  seurat.rna = subset(seurat.rna, projCtype == stage_type)
  
  mtx.mllr = seurat.rna[['RNA']]@counts
  
  efreq.age.mllr <- sapply(c('Young', 'Old'), function(x) {
    
    cl_data <- mtx.mllr[, seurat.rna$Age_Group == x]
    
    Matrix::rowMeans(cl_data > 0)
    
  })
  
  save(efreq.age.mllr, file = paste0('MetaData/scRNA/efreq_ctype_YoungvsOld_',
                                     stage_type, '.RData'))
  
  ann.mllr = ann.mllr[Ctype0 == 'Blasts']
  cls = subset(ann.mllr, select = c(Age_Group, sample_cell))
  cls = subset(ann.mllr, projCtype == stage_type,
                         select = c(Age_Group, sample_cell))
  names(cls) = c('cluster', 'barcode')
  dscore.mllr[is.na(dscore.mllr)] = 0
  
  cls = cls[barcode %in% cnames]
  mcell = 2000
  if(stage_type %in% c('Early-Prog-like', 'Immature-B-like', 'Mature-B-like'))
    mcell = 100
  
  tf.mllr.sele <- runDiffMotifEnrich(dscore.mllr, clusters = cls, fdr = 0.05,
                                max_cell_per_clust = mcell,
                                topn = 200, min_frac_per_cluster = 0.2)
  
  sele.tfs0 = lapply(c('Young', 'Old'), function(x){
    res.diff0 = tf.mllr.sele[order(-mean1)]
    tfs = res.diff0[cluster1 == x]$feature
    efreq0 = efreq.age.mllr[, x]
    ## filter by expr frac
    tfs = tfs[tfs %in% names(efreq0)]
    tfs = tfs[efreq0[tfs] > expr_cutoff]
    
    return(res.diff0[feature %in% tfs, ])
  })
  names(sele.tfs0) = c('Young', 'Old')
  if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('rbind', sele.tfs0)
  
  tf.mllr.sele = tf.mllr.sele[feature %in% sele.tfs0$feature]
  tf.mllr.sele = tf.mllr.sele[, delta := (mean1 - mean0)]
  tf.mllr.sele = tf.mllr.sele[order(-delta), ]
 
  nshow = min(20, nrow(tf.mllr.sele))  ## update with each group 10 (previous totally20)
  
  tf.mllr.sele = tf.mllr.sele[1:nshow, ]
  tf.res.list[[stage_type]] = tf.mllr.sele
  
  age_color <- c( 
    'Young' = '#A6CEE3',
    'Old' = '#1F78B4'
  )
  ndownsample = min(1000, nrow(cls))
  
  ## plot group mean
  zscore.mllr1 = zscore.mllr[, cls$barcode]
  zscore.mllr1 = t(scale(t(zscore.mllr), center = T, scale = T))
  
  avg.y <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, 
                                 cls[cluster=='Young']$barcode])
  avg.o <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, 
                                 cls[cluster=='Old']$barcode])
  
  avg.score = cbind(avg.y, avg.o)
  colnames(avg.score) = c('Young', 'Old')
  avg.score.list[[stage_type]] = avg.score
  #pheatmap::pheatmap(avg.score, cluster_cols = F, 
  #                                            cluster_rows = T)
}
saveRDS(tf.res.list, file = paste0('EP_Prediction/Blasts_YoungvsOld_enriched_tfs_list_by_stage_exprfrac', 
                                   expr_cutoff, '.rds'))

tfs = lapply(avg.score.list, function(x) rownames(x) )
tfs = unique(do.call('c', tfs))
ann.mllr0 = ann.mllr[sample_cell %in% colnames(zscore.mllr)]
zscore.mllr0 = zscore.mllr[, ann.mllr0$sample_cell]
zscore.mllr0 = zscore.mllr0[tfs, ]

dscore.mllr0 = dscore.mllr[, ann.mllr0$sample_cell]
dscore.mllr0 = dscore.mllr0[tfs, ]


stages = c('Early-Prog-like', 'Pre-pro-B-like', 'Pro-B-like',
           'Pre-B-like', 'Immature-B-like', 'Mature-B-like')
zscore.diff <- sapply(stages, function(x){
  bc1 = ann.mllr0[projCtype == x & Age_Group == 'Young']$sample_cell
  bc0 = ann.mllr0[projCtype == x & Age_Group == 'Old']$sample_cell
  
  z1 <- rowMeans(zscore.mllr0[, bc1])
  z0 <- rowMeans(zscore.mllr0[, bc0])
  return(z1 - z0)
})

dscore.diff <- sapply(stages, function(x){
  bc1 = ann.mllr0[projCtype == x & Age_Group == 'Young']$sample_cell
  bc0 = ann.mllr0[projCtype == x & Age_Group == 'Old']$sample_cell
  
  d1 <- rowMeans(dscore.mllr0[, bc1])
  d0 <- rowMeans(dscore.mllr0[, bc0])
  return(d1 - d0)
})

ann_column = data.frame('stage' = stages,
                        'barcode' = stages,
                        stringsAsFactors = F)
rownames(ann_column) = stages
ann_column$barcode <- NULL
ann_color = list('stage' = c("Early-Prog-like" = "#E3191C", "Pre-pro-B-like" = "#347EB5",
              "Pro-B-like" = "#4CAC47", "Pre-B-like" = "#A35728", 
              "Immature-B-like" = "#924C9F", "Mature-B-like" = "#00FFFF"))
zscore.tailor <- zscore.diff

ph_bulk <- pheatmap::pheatmap(zscore.tailor, cluster_cols = F, 
                   cluster_rows = T, annotation_col = ann_column,
                   annotation_colors = ann_color, show_colnames = F,
                   border_color = "NA", 
                   breaks = sort(unique(c(seq(-1, 0, length.out = 50), 
                              seq(0, 5, length.out = 50)))),
                   color = colorRampPalette(c("blue", "white", "orange"))(99))

hmfile_avg = paste0(output_dir, '/enriched_TFs_YoungVsOld_AllStages_zscoreAvgDiff.eps')
ggsave(ph_bulk, filename = hmfile_avg,
       width = 8, height = 8, device = 'eps')

## put the non-significant one's as zero (in the heatmap)
tf_stage_map = matrix(0, length(tfs), length(stages))
rownames(tf_stage_map) = tfs
colnames(tf_stage_map) = stages
for(stage0 in stages){
  ffs = unique(tf.res.list[[stage0]]$feature)
  tf_stage_map[ffs, stage0] = 1
}

ph_bulk <- pheatmap::pheatmap(zscore.tailor * tf_stage_map, cluster_cols = F, 
                              cluster_rows = T, annotation_col = ann_column,
                              annotation_colors = ann_color, show_colnames = F,
                              border_color = "NA", 
                              breaks = sort(unique(c(seq(-1, 0, length.out = 50), 
                                                     seq(0, 4, length.out = 50)))),
                              color = colorRampPalette(c("blue", "white", "orange"))(99))
hmfile_avg = paste0(output_dir, '/enriched_TFs_YoungVsOld_AllStages_zscoreAvgDiff_updated.eps')
ggsave(ph_bulk, filename = hmfile_avg,
       width = 8, height = 8, device = 'eps')


## write into table as supple table 
for(stage0 in names(tf.res.list)){
  dd0 = tf.res.list[[stage0]]
  dd0[, 'diff_dev' := ifelse(cluster1 == 'Old', 
                             mean0 - mean1, mean1 - mean0 )]
  dd0 = subset(dd0, select=c("feature", "diff_dev", "pv_adjust"))
  names(dd0) = c('TF', 'diff_dev(Young-old)', 'FDR')
  write.table(dd0, file = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac0.1/Tables/tables_enriched_tf_', 
                                 stage0, '.txt'), row.names = F, sep = '\t', quote = F)
}


## update B- vs M-lineage comparison ####
## prepare dscore and zscore
load('chromVAR_Objects/zscore_dscore_18MLLr_allPeaks.RData')
ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
ann.mllr = ann.mllr[Ctype0 == 'Blasts']
seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna.mllr = subset(seurat.rna.mllr, Ctype0 == 'Blasts')
seurat.rna.mllr$lineage = 'Others'

seurat.rna.mllr$lineage[seurat.rna.mllr$projCtype %in% c('CLP-like', 'Pre-pro-B-like',
                                                 'Pro-B-like', 'Pre-B-like',
                                                 'Immature-B-like', 'Mature-B-like')] = 'B_lineage'
seurat.rna.mllr$lineage[seurat.rna.mllr$projCtype %in% c('cDC-like', 'DC-Progenitor-like',
                                                 'pDC-like',
                                                 'GMP-like', 'Mono-like')] = 'M_lineage'
efreq.lineage <- sapply(c('B_lineage', 'M_lineage'), function(x){
  return(rowMeans(seurat.rna.mllr@assays$RNA@counts[, seurat.rna.mllr$lineage == x] > 0))
})

rm(seurat.rna.mllr)
save(efreq.lineage, file = 'MetaData/scRNA/efreq_ctype_BvsM.RData')

load('MetaData/scRNA/efreq_ctype_BvsM.RData')
ann.mllr$lineage = 'Others'
ann.mllr$lineage[ann.mllr$projCtype %in% c('CLP-like', 'Pre-pro-B-like',
                                                 'Pro-B-like', 'Pre-B-like',
                                                 'Immature-B-like', 'Mature-B-like')] = 'B_lineage'
ann.mllr$lineage[ann.mllr$projCtype %in% c('cDC-like', 'DC-Progenitor-like',
                                                 'pDC-like',
                                                 'GMP-like', 'Mono-like')] = 'M_lineage'


## change cell names to be consistent
cls = subset(ann.mllr, select = c(sample_cell, lineage))
cls = data.table(cls, keep.rownames = T)
names(cls)[1:2] = c('barcode', 'cluster')
cls = cls[cluster != 'Others']

cnames = colnames(dscore.mllr)
cnames = sapply(cnames, function(x) gsub('_scATAC', '', x))
cnames = sapply(cnames, function(x) gsub('MLL_', 'MLLr', x))
dscore.mllr = dscore.mllr[, cnames %in% cls$barcode]
zscore.mllr = zscore.mllr[, cnames %in% cls$barcode]
colnames(dscore.mllr) = cnames[cnames %in% cls$barcode]
colnames(zscore.mllr) = cnames[cnames %in% cls$barcode]

expr_cutoff = 0.1 
output_dir = paste0('Figures/scATAC/TF_Enrich/MLLr18Infants/filteredByExpr_frac', 
                    expr_cutoff)

dscore.mllr[is.na(dscore.mllr)] = 0


cls = cls[barcode %in% cnames]
mcell = 2000
tf.mllr <- runDiffMotifEnrich(dscore.mllr, clusters = cls, fdr = 0.05,
                              max_cell_per_clust = mcell,
                              topn = 200, min_frac_per_cluster = 0.2)

sele.tfs0 = lapply(c('B_lineage', 'M_lineage'), function(x){
  res.diff0 = tf.mllr[order(-mean1)]
  tfs = res.diff0[cluster1 == x]$feature
  efreq0 = efreq.lineage[, x]
  ## filter by expr frac
  tfs = tfs[tfs %in% names(efreq0)]
  tfs = tfs[efreq0[tfs] > expr_cutoff]
  
  return(res.diff0[feature %in% tfs, ])
})
names(sele.tfs0) = c('B_lineage', 'M_lineage')
saveRDS(sele.tfs0, "EP_Prediction/Blasts_BvsM_enriched_tfs_exprfrac0.1.rds")

if(class(sele.tfs0) == 'list') sele.tfs0 = do.call('rbind', sele.tfs0)


tf.mllr.sele = tf.mllr[feature %in% sele.tfs0$feature]
tf.mllr.sele = tf.mllr.sele[, delta := (mean1 - mean0)]
tf.mllr.sele = tf.mllr.sele[order(-delta), ]

tf1 = tf.mllr.sele[cluster1 == 'M_lineage', ]
tf2 = tf.mllr.sele[cluster1 == 'B_lineage', ]

tf.mllr.sele = rbind(tf1[1:20, ], tf2[1:10, ])

age_color <- c( 
  'B_lineage' = '#A6CEE3',
  'M_lineage' = '#1F78B4'
)
ndownsample = min(1000, nrow(cls))

ph <- plot_enrich_tf1(unique(tf.mllr.sele$feature), zscore.mllr,
                      bc_clusters = cls[cls$cluster != 'Others', ],
                      cluster.levels =c('M_lineage', 'B_lineage'), 
                      cluster.col = age_color,
                      up.qt = 0.94, low.qt = 0.02,
                      ndownsample =  ndownsample, match_cell = T,
                      reScale = F, cluster.rows = F,
                      color_style = 'purple-yellow')

hmfile = paste0(output_dir, '/enriched_TFs_BVsM.eps')
ht = ceiling(12 * nrow(tf.mllr.sele)/ 50 )
ggsave(ph, filename = hmfile,
       width = 12, height = 12, device = 'eps')

## plot group mean
zscore.mllr1 = zscore.mllr[, cls$barcode]
zscore.mllr1 = t(scale(t(zscore.mllr), center = T, scale = T))

avg.y <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, cls[cluster=='B_lineage']$barcode])
avg.o <- rowMeans(zscore.mllr1[tf.mllr.sele$feature, cls[cluster=='M_lineage']$barcode])

avg.score = cbind(avg.y, avg.o)
colnames(avg.score) = c('B_lineage', 'M_lineage')
avg.score[avg.score > 1.5] = 1.5
ph_bulk <- pheatmap::pheatmap(avg.score, cluster_cols = F, 
                              cluster_rows = T)

hmfile_avg = paste0(output_dir, '/enriched_TFs_BvsM_avg.eps')
ggsave(ph_bulk, filename = hmfile_avg,
       width = 4, height = ht, device = 'eps')



## plot example for HSPC, HSPC1 and Blasts ####
load('chromVAR_Objects/zscore_dscore_all.RData')
ann.hd = fread('MetaData/scATAC/hd_bc_sample_inf.txt')
ann.mllr = fread('MetaData/scATAC/mllr18Infants_bc_sample_inf.txt')
ann.hspc1 = ann.mllr[Ctype0 == 'Progenitors']
ann.blast = ann.mllr[Ctype0 == 'Blasts']
ann.hspc = ann.hd[Ctype == 'HSPC']

dscore.sele = dscore.all[, c(ann.hspc$sample_cell, ann.hspc1$sample_cell,
                              ann.blast$sample_cell)]
dt = data.frame('Ctype' = rep(c('HSPC', 'HSPC1', 'Blasts'),
                              c(nrow(ann.hspc), nrow(ann.hspc1), nrow(ann.blast))),
                stringsAsFactors = F)

rownames(dt) = colnames(dscore.sele)

seurat.sele = CreateSeuratObject(dscore.sele, assay = 'TF',
                                 meta.data = dt)
seurat.sele$Ctype = factor(seurat.sele$Ctype, levels = c('HSPC', 'HSPC1', 'Blasts'))
myColor = c("#A6CEE3", "#1F78B4", "red")

p0 <- StackedVlnPlot(seurat.sele, features = c('GATA2', 'HMGA2', 'RUNX1', 'HOXA9',
                                         'MEIS1', 'MEF2C', 'EBF1', 'PAX5'),
               group.by = 'Ctype', myColors = myColor) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(p0, filename = 'Figures/summary/example_tfs_vlnplot.eps', height = 6,
       width = 3, device = 'eps')

dt$bc = rownames(dt)
mtx <- sapply(c('HSPC', 'HSPC1', 'Blasts'), function(x){
  rowMeans(dscore.sele[, rownames(dt[dt$Ctype == x, ])])
})

p1 <- pheatmap::pheatmap(mtx[c('GATA2', 'HMGA2', 'RUNX1', 'HOXA9',
                       'MEIS1', 'MEF2C', 'EBF1', 'PAX5'), ],
                   cluster_cols = F, 
                   cluster_rows = T)
ggsave(p1, filename = 'Figures/summary/example_tfs_hm.eps', height = 6,
       width = 4, device = 'eps')
