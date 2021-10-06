source('scDataAnalysis_Utilities.R')

##************************************************##
## -- sum cell-cell-communication results ##
## -- using overlap of 2/3 methods ##
##************************************************##

## list of functions

dot_plot.cellphonedb = function(selected_rows = NULL,
                                selected_columns = NULL,
                                filename = 'plot.pdf',
                                width = 8,
                                height = 10,
                                means_path = './means.txt',
                                pvalues_path = './pvalues.txt',
                                means_separator = '\t',
                                pvalues_separator = '\t',
                                output_extension = '.eps', 
                                save_plot = T){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, 
                        sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, 
                         sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = min(pr[pr!=0])/10
  plot.data = cbind(plot.data, 1+log10(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  plot.data = plot.data[order(factor(plot.data$clusters,
                                     levels = selected_columns)), ] 
  plot.data$clusters = sapply(plot.data$clusters, function(x) gsub('Monocyte', 'Mono', x))
  plot.data$clusters = sapply(plot.data$clusters, function(x) gsub('NK_T', 'NKT', x))
  
  selected_columns = sapply(selected_columns, function(x) gsub('Monocyte', 'Mono', x))
  selected_columns = sapply(selected_columns, function(x) gsub('NK_T', 'NKT', x))
  
  
  plot.data$clusters = factor(plot.data$clusters, levels = selected_columns)
  
  ## only keep significant pathways for given clusters
  plot.data = data.table(plot.data)
  plot.data[, 'minp' := min(pvalue), by = pair]
  plot.data = plot.data[minp <= 0.05]
  
  #  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  my_palette <- viridis(100)
  pdot <- ggplot(plot.data, aes(x=clusters, y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log10 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, family = 'Arial'),
          axis.text.y = element_text(size=16, colour = "black", family = 'Arial'),
          axis.title=element_blank(),
          text = element_text('Arial'),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black")) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
  
  if (output_extension == '.esp') {
    if(save_plot) ggsave(plot = pdot, filename, width = width, height = height,
                         device = 'eps', limitsize=F)
  }
  else {
    if(save_plot) ggsave(plot = pdot, filename, width = width, height = height, 
                         limitsize=F)
  }
  plot.data$pair = as.character(plot.data$pair)
  return(list('plot' = pdot, 'pathways' = unique(plot.data$pair), 
              'plot.data' = plot.data))
}

## plot interactions for an interested ctype
## interestCtype: centered ctype
## otherCtypeGroup: normal/blasts/all
## interestTerms: NULL--plot all sig terms
## ctype.pair.order: left: interectCtype | other; right: other|interestCtype,
## two-side: both left and right
dotPlot4SingleCtype.cellphondb <- function(sampleName, interestCtype = 'HSPC1', 
                                           otherCtypes = norm.types,
                                           interestTerms.key = 'IFNG', 
                                           plotName.key = 'HSPC1vsImmnue',
                                           ctype.pair.order = 'two-side',
                                           interestTerms.order = 'right',
                                           save_plot1 = T, save_plot2 = T){
  cb_outdir = paste0('CellPhoneDB_18MLLr/data4Yuxuan/', sampleName)
  
  if(any(interestCtype == 'Blasts') || any(otherCtypes == 'Blasts'))
    cb_outdir = paste0('CellPhoneDB_18MLLr/data4Yuxuan/', sampleName, '/Blasts')
  
  dd = fread(paste0(cb_outdir, '/significant_means.txt'))
  #fig_outdir = 'CellPhoneDB_18MLLr/Figures'
  fig_outdir = cb_outdir
  nn = names(dd)
  nn0 = nn[grepl(nn, pattern = interestCtype)]
  if(length(nn0) == 0) return(NULL)
  dd0 = subset(dd, select = nn0)
  dd0 = sapply(1:nrow(dd0), function(x) if(!all(is.na(dd0[x, ]))) return(dd0[x, ]) )
  len0 = sapply(dd0, length)
  
  dd.sele = dd[which(len0 > 0), ]
  
  
  dd.sele = subset(dd.sele, select = c(nn[1:12], nn0))
  hspc.pair = subset(dd.sele, select = nn0)
  dd.sele$ctype_pair <- sapply(1:nrow(hspc.pair), function(x) 
    paste(nn0[which(!is.na(hspc.pair[x,]))], collapse = ','))
  
  dd.sele = subset(dd.sele, select = c(nn[1:12], 'ctype_pair'))
  
  dd.sele[, 'gene_pair' := paste0(gene_a, ';', gene_b)]
  
  dd.sele = dd.sele[order(rank)]
  sele.rows = dd.sele$interacting_pair
  
  ctype.pairs <- data.table('ctype1' = interestCtype, 
                            ctype2 = setdiff(otherCtypes, interestCtype))
  left.pairs = paste0(ctype.pairs$ctype1, '|', ctype.pairs$ctype2)
  right.pairs = paste0(ctype.pairs$ctype2, '|', ctype.pairs$ctype1)
  
  if(ctype.pair.order == 'two-side'){
    #ctype.pairs[, 'cpair' := paste0(paste0(ctype1, '|', ctype2), ',', 
    #                              paste0(ctype2, '|', ctype1))]
    #sele.cols = c(sapply(ctype.pairs$cpair, function(x) unlist(strsplit(x, ','))), 
    #             paste0(interestCtype, '|', interestCtype))
    
    sele.cols = c(left.pairs, right.pairs)
    sele.cols = sele.cols[sele.cols%in%nn0]
  }
  if(ctype.pair.order == 'left'){
    sele.cols = left.pairs[left.pairs%in%nn0]
  }
  if(ctype.pair.order == 'right'){
    sele.cols = right.pairs
    sele.cols = right.pairs[right.pairs%in%nn0]
  }
  
  ## which type of columns are plotted
  filename = paste0(fig_outdir,  '/', sampleName, '_dotplot_interaction_relate2', 
                    plotName.key, '.eps')
  
  if(!is.null(interestTerms.key)) filename0 = paste0(fig_outdir, '/', sampleName, '_dotplot_interaction_relate2', 
                                                     plotName.key, '_4', interestTerms.key, '.eps')
  
  
  res <- dot_plot.cellphonedb(selected_rows = sele.rows, selected_columns = sele.cols,
                              means_path = paste0(cb_outdir, '/means.txt'),
                              pvalues_path = paste0(cb_outdir, '/pvalues.txt'),
                              filename = filename, save_plot = save_plot1,
                              width = 16, height = 18, output_extension = '.eps')
  
  if(!is.null(interestTerms.key)) {
    
    interest.terms = sele.rows[grepl(sele.rows, pattern = interestTerms.key)]
    
    # only plot interested terms
    if(length(interest.terms) > 0) {
      if(interestTerms.order == 'left') sele.cols = left.pairs[left.pairs%in%nn0]
      if(interestTerms.order == 'right') sele.cols = right.pairs[right.pairs%in%nn0]
      
      res = dot_plot.cellphonedb(selected_rows = interest.terms, selected_columns = sele.cols,
                                 means_path = paste0(cb_outdir, '/means.txt'),
                                 pvalues_path = paste0(cb_outdir, '/pvalues.txt'),
                                 filename = filename0, 
                                 width = 16, height = 4, output_extension = '.eps', 
                                 save_plot = save_plot2)
      
    }
  }
  return(res)
}

## adjust some names to enable overlapping
adjust_cb_intName <- function(pdata){
  ## IFNG
  pdata0 = pdata[grepl(pair, pattern = 'IFNG')]
  pdata1 = pdata[!grepl(pair, pattern = 'IFNG')]
  if(nrow(pdata0) > 0){
    tmp1 = tmp2 = pdata0
    tmp1$pair = 'IFNG_IFNGR1'
    tmp2$pair = 'IFNG_IFNGR2'
    pdata = rbind(pdata1, tmp1, tmp2)
  }
  
  ## 
  pdata[, 'pair' := gsub('TGFbeta receptor', 'TGFBR', pair, fixed = T)]
  pdata[, 'pair' := gsub(' receptor', 'R', pair, fixed = T)]
  
  ## aLb2 complex
  pdata0 = pdata[grepl(pair, pattern = 'aLb2')]
  pdata1 = pdata[!grepl(pair, pattern = 'aLb2')]
  if(nrow(pdata0) > 0){
    tmp1 = tmp2 = pdata0
    tmp1$pair = sapply(tmp1$pair, function(x) gsub('aLb2 complex', 'ITGAL', x))
    tmp2$pair = sapply(tmp2$pair, function(x) gsub('aLb2 complex', 'ITGB2', x))
    pdata = rbind(pdata1, tmp1, tmp2)
  }
  
  ## a4b1 complex
  pdata0 = pdata[grepl(pair, pattern = 'a4b1')]
  pdata1 = pdata[!grepl(pair, pattern = 'a4b1')]
  if(nrow(pdata0) > 0){
    tmp1 = tmp2 = pdata0
    tmp1$pair = sapply(tmp1$pair, function(x) gsub('a4b1 complex', 'ITGB1', x))
    tmp2$pair = sapply(tmp2$pair, function(x) gsub('a4b1 complex', 'ITGA4', x))
    pdata = rbind(pdata1, tmp1, tmp2)
  }
  return(pdata)
}

## annotate active/repressive interactions ####
ann.ptw = c('CD160_TNFRSF14' = 1, 'KLRC1_HLA-E' = -1,
            'TGFBR3_TGFB1' = -1, 'TGFBR2_TGFB1' = -1,
            'TGFBR1_TGFB1' = -1, 'TGFB1_TGFBR1' = -1, 
            'SELPLG_SELL' = -1,
            'SELL_SELPLG' = -1, 'TNFRSF1B_TNF' = -1,
            'SELL_CD34' = -1, 'HAVCR2_LGALS9' = -1,
            'ADRB2_IL1B' = -1, 'KIR2DL3_HLA-C' = -1,
            'LTB_LTBR' = -1, 'CD226_NECTIN2' = 1,
            'CD44_HGF' = -1, 'TIGIT_NECTIN2' = -1,
            'KLRC2_HLA-E' = 1, 'CD2_CD58' = 1,
            'CD58_CD2' = 1,
            'CD244_CD48' = 1, 'CD48_CD244' = 1,
            'ITGAL_ICAM2' = 1,
            'ITGAL_ICAM3' = 1, 'ITGAL_ICAM4' = 1,
            'ITGB2_ICAM2' = 1, 'ITGB2_ICAM3' = 1,
            'ITGB2_ICAM4' = 1, 'CCR5_CCL5' = -1,
            'CD44_HBEGF' = -1, 'TNFRSF1A_TNF' = -1,
            "TNFRSF14_TNFSF14" = -1, "TNFSF14_TNFRSF14" = -1,
            'ICAM2_ITGB2' = 1, 'ICAM3_ITGAL' = 1,
            'ICAM3_ITGB2' = 1, 'ICAM2_ITGAL' = 1,
            'SEMA4D_PLXNB2'= -1, 'SEMA4D_CD72' = 1,
            'FLT3_FLT3LG' = 1, 'FLT3LG_FLT3' = -1,
            'IFNG_IFNGR1' = 1, 'IFNG_IFNGR2' = 1)

ann.ptw.nk = c('CD160_TNFRSF14' = 1, 'KLRC1_HLA-E' = -1,
               'TGFBR3_TGFB1' = -1, 'TGFBR2_TGFB1' = -1,
               'TGFBR1_TGFB1' = -1, 'TGFB1_TGFBR1' = -1, 
               'SELPLG_SELL' = 0,
               'SELL_SELPLG' = 0, 'TNFRSF1B_TNF' = -1,
               'SELL_CD34' = 0, 'HAVCR2_LGALS9' = -1,
               'ADRB2_IL1B' = 0, 'KIR2DL3_HLA-C' = -1,
               'LTB_LTBR' = -1, 'CD226_NECTIN2' = 1,
               'CD44_HGF' = -1, 'TIGIT_NECTIN2' = -1,
               'KLRC2_HLA-E' = 1, 'CD2_CD58' = 1,
               'CD58_CD2' = 1,
               'CD244_CD48' = 1, 'CD48_CD244' = 1,
               'ITGAL_ICAM2' = 1,
               'ITGAL_ICAM3' = 1, 'ITGAL_ICAM4' = 1,
               'ITGB2_ICAM2' = 1, 'ITGB2_ICAM3' = 1,
               'ITGB2_ICAM4' = 1, 'CCR5_CCL5' = -1,
               'CD44_HBEGF' = -1, 'TNFRSF1A_TNF' = -1,
               "TNFRSF14_TNFSF14" = -1, "TNFSF14_TNFRSF14" = -1,
               'ICAM2_ITGB2' = 1, 'ICAM3_ITGAL' = 1,
               'ICAM3_ITGB2' = 1, 'ICAM2_ITGAL' = 1,
               'SEMA4D_PLXNB2'= 0, 'SEMA4D_CD72' = 1,
               'FLT3_FLT3LG' = -1, 'FLT3LG_FLT3' = -1,
               'IFNG_IFNGR1' = 1, 'IFNG_IFNGR2' = 1)

## prepare DEGs between HSPC1 and Blasts ####
seurat.rna.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna.mllr = subset(seurat.rna.mllr, Ctype0 %in% c('Blasts', 'Progenitors'))
Idents(seurat.rna.mllr) = seurat.rna.mllr$Ctype0
degs = FindAllMarkers(seurat.rna.mllr, only.pos = T, 
                      max.cells.per.ident = 500,
                      logfc.threshold = 0.02)
degs = data.table(degs)
write.table(degs, file = 'MetaData/scRNA/degs_HSPC1vsBlasts.txt',
            row.names = F, sep = '\t', quote = F)


## < sum for HSPC1 ####
degs = fread('MetaData/scRNA/degs_HSPC1vsBlasts.txt')
topn = 0.15
seleImmune.types = c('CD8_Cytotoxic', 'CD16_NK', 'CD56_NK')
others.ctype.name = 'selectedImmune'
others.ctype = switch (others.ctype.name,
                       'selectedImmune' = seleImmune.types,
                       'allImmune' = norm.types)

data.name = 'PoolPatientsWHSPC1'

interested.ctype = 'HSPC1'
res1 = dotPlot4SingleCtype.cellphondb(sampleName = data.name, 
                                      interestCtype = interested.ctype, 
                                      otherCtypes = others.ctype, 
                                      interestTerms.key = NULL,
                                      plotName.key = paste0(interested.ctype, 'vsImmune'),
                                      ctype.pair.order = 'left')
res2 = dotPlot4SingleCtype.cellphondb(sampleName = data.name, 
                                      interestCtype = interested.ctype, 
                                      otherCtypes = others.ctype, 
                                      interestTerms.key = NULL,
                                      plotName.key = paste0('ImmnuneVS', interested.ctype),
                                      ctype.pair.order = 'right')

pdata2 = subset(res2$plot.data, select = c(pair, clusters, pvalue, mean))
pdata1 = subset(res1$plot.data, select = c(pair, clusters, pvalue, mean))

## rename interaction and ctype pair to the order of immune vs HSPC1
pdata1$pair = paste0(sapply(pdata1$pair, function(x) unlist(strsplit(x, '_'))[2]), '_',
                     sapply(pdata1$pair, function(x) unlist(strsplit(x, '_'))[1]))
pdata1$clusters = as.character(pdata1$clusters)
pdata2$clusters = as.character(pdata2$clusters)
pdata1$clusters = paste0(sapply(pdata1$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[2]),
                         '|',
                         sapply(pdata1$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[1]))

pdata = rbind(pdata1, pdata2)
pdata = pdata[pvalue <= 0.01]



pdata = adjust_cb_intName(pdata)
#pdata[, 'clusters' := gsub('Mono|', 'Monocyte|', clusters, fixed = T)]

## filter by KumarMethod
dd.k <- fread(paste0('Cell_Interactions_Sum/Kumar/', data.name, '/interaction_', interested.ctype,
                     '_Immune_top', topn*100, 'Perc.txt'))
pdata.k = NULL
for(cluster0 in unique(pdata$clusters)){
  pdata0 = pdata[clusters == cluster0]
  dd.k0 = dd.k[clusters == cluster0]
  pdata0 = pdata0[pair %in% dd.k0$pair]
  pdata.k = rbind(pdata.k, pdata0)
}


## further filter by DEGs 
names(degs)[which(names(degs) == 'gene')] = 'gene_name'
degs = degs[p_val < 0.05]

pdata.k$gene1 = sapply(pdata.k$pair, function(x) unlist(strsplit(x, '_'))[1])
pdata.k$gene2 = sapply(pdata.k$pair, function(x) unlist(strsplit(x, '_'))[2])
pdata.k = pdata.k[gene2 %in% degs[cluster == 'Progenitors']$gene_name]

pdata.k[clusters == 'CD8_Cytotoxic|HSPC1', 
        'direction' := ifelse(pair %in% names(ann.ptw), ann.ptw[pair], 0)]
pdata.k[clusters != 'CD8_Cytotoxic|HSPC1', 
        'direction' := ifelse(pair %in% names(ann.ptw.nk), ann.ptw.nk[pair], 0)]


## plot overlapped results with kumarMethod
pdata.k$isSig = 1
pdata.mtx = reshape2::recast(pdata.k, pair ~ clusters, measure.var = 'isSig')
rownames(pdata.mtx) = pdata.mtx[, 1]
pdata.mtx = pdata.mtx[, -1]
pdata.mtx[is.na(pdata.mtx)] = 0
pdata.mtx = as.matrix(pdata.mtx)

#pdata.mtx = -log10(pdata.mtx[, -1])
#pdata.mtx[pdata.mtx <= 1] = 1
ph <- pheatmap::pheatmap(pdata.mtx, cluster_cols = T, cluster_rows = T,
                         angle_col = 45, scale = 'none', 
                         clustering_distance_rows = 'euclidean',
                         color = colorRampPalette(c( "white", "#DC143C"))(100))
fname = paste0('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/', data.name, '_relate2_',
               interested.ctype, '_vs_', others.ctype.name, '_heatmap.eps')
ggsave(ph, filename = fname,
       width = 14, height = 7, device = 'eps')


pdata.k = pdata.k[order(-direction)]
pdata.k$pair = factor(pdata.k$pair, levels = unique(pdata.k$pair))
pdata.k$clusters = factor(pdata.k$clusters, levels = c("CD8_Cytotoxic|HSPC1",
                                                       "CD16_NK|HSPC1", "CD56_NK|HSPC1"))

names(pdata.k)[which(names(pdata.k) == 'mean')] = 'expression'

pd.cyto <- ggplot(pdata.k, aes(x=clusters, y=pair)) +
  geom_point(aes(size=expression, color=expression)) +
  scale_color_gradientn('expression', 
                        colors = viridis(500)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, family = 'Arial'),
        axis.text.y = element_text(size=14, colour = "black", family = 'Arial',
                                   margin = margin(t = 0, r = 0, b = 0, l = 20)),
        axis.title=element_blank(),
        text = element_text('Arial'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 

fname.dot1 = paste0('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/', data.name, '_relate2_',
                    interested.ctype, '_vs_', others.ctype.name, '_dotplot_cytotoxic.eps')
ggsave(pd.cyto, filename = fname.dot1,
       width = 7, height = 7, device = 'eps')


## trying circlize plot
if(T){
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  myColors = getPalette(15)
  names(myColors) = c('CD14_Monocyte', 'CD16_Monocyte', 'CD16_NK', 'CD4_Early_Activation', 'CD4_Naive', 'CD56_NK', 'CD8_Cytotoxic', 'CD8_EffectorMemory',
                      'CD8_Naive', 'cDC', 'Exhausted_T', 'Mature_B', 'NKT', 'pDC', 'Unknown')
  
  pdata.k$clusters = as.character(pdata.k$clusters)
  pdata.k$ctype1 = sapply(pdata.k$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[1])
  pdata.k$ctype2 = sapply(pdata.k$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[2])
  
  pdata.c <- subset(pdata.k, select = c(gene1, gene2, ctype1, ctype2,
                                                                        expression, expression, direction))
  names(pdata.c) = c("receptor", "ligand", "cell_to", "cell_from" , 
                     "cell_to_mean_exprs", "cell_from_mean_exprs", 'direction')
  pdata.c$comm_type = 'checkpoint'
  
  
  cell.col = myColors[unique(pdata.c$cell_to)]
  cell.col['HSPC1'] = "#984EA3"
  
  pdata.c$link.col = 'gray70'
  pdata.c$link.col = ifelse(pdata.c$direction == -1, 'red',
                            pdata.c$link.col)
  pdata.c$link.col = ifelse(pdata.c$direction == 1, 'skyblue',
                            pdata.c$link.col)
  pdata.c = pdata.c[order(direction)]
  postscript('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/circular_plot_HSPC1vsImmune.eps',
             height = 6, width = 6)
  iTALK::LRPlot(pdata.c, datatype = 'mean count', 
                cell_col = cell.col, transparency = 0,
                link.arr.width = 0.2,
                link.arr.lwd = pdata.c$cell_to_mean_exprs,
                link.arr.col = pdata.c$link.col)
  dev.off()
  
}


## killing score
if(T){
  pdata.k[, 'kscore' := sum(expression * direction), by = clusters]
  kscore <- pdata.k %>% subset(select = c(clusters, kscore)) %>% .[!duplicated(.)]
  
  kcolor = c('CD8_Cytotoxic|HSPC1' = "#C4625D",
             'CD16_NK|HSPC1' = "#3A85A8",
             'CD56_NK|HSPC1' = "#8D5B96")
  p0 <- ggplot(data = kscore, 
               aes(x = clusters, y = kscore, fill = clusters)) +
    geom_bar(stat = 'identity') + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    NoLegend() + xlab('') + ylab('Pro-killing score') +
    scale_fill_manual(values = kcolor)
  
  write.table(kscore, 
              file = 'Figures/summary/Cell_cell_communication/HSPC1vsBlasts/cytotoxic_pro_killing_score_HSPC1.txt',
              row.names = F, quote = F, sep = '\t')
}

## < sum for corresponding Blasts ####
degs = fread('MetaData/scRNA/degs_HSPC1vsBlasts.txt')
blast.types = c('Blasts')
seleImmune.types = c('CD8_Cytotoxic', 'CD16_NK', 'CD56_NK')

others.ctype.name = 'selectedImmune'

others.ctype = switch (others.ctype.name,
                       'selectedImmune' = seleImmune.types,
                       'allImmune' = norm.types)

data.name = 'PoolPatientsWHSPC1'

interested.ctype = 'Blasts'
res1 = dotPlot4SingleCtype.cellphondb(sampleName = data.name, 
                                      interestCtype = interested.ctype, 
                                      otherCtypes = others.ctype, 
                                      interestTerms.key = NULL,
                                      plotName.key = paste0(interested.ctype, 'vsImmune'),
                                      ctype.pair.order = 'left')
res2 = dotPlot4SingleCtype.cellphondb(sampleName = data.name, 
                                      interestCtype = interested.ctype, 
                                      otherCtypes = others.ctype, 
                                      interestTerms.key = NULL,
                                      plotName.key = paste0('ImmnuneVS', interested.ctype),
                                      ctype.pair.order = 'right')

pdata2 = subset(res2$plot.data, select = c(pair, clusters, pvalue, mean))
pdata1 = subset(res1$plot.data, select = c(pair, clusters, pvalue, mean))

## rename interaction and ctype pair to the order of immune vs HSPC1
pdata1$pair = paste0(sapply(pdata1$pair, function(x) unlist(strsplit(x, '_'))[2]), '_',
                     sapply(pdata1$pair, function(x) unlist(strsplit(x, '_'))[1]))
pdata1$clusters = as.character(pdata1$clusters)
pdata2$clusters = as.character(pdata2$clusters)
pdata1$clusters = paste0(sapply(pdata1$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[2]),
                         '|',
                         sapply(pdata1$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[1]))

pdata = rbind(pdata1, pdata2)
pdata = pdata[pvalue <= 0.01]


pdata = adjust_cb_intName(pdata)
pdata[, 'clusters' := gsub('Mono|', 'Monocyte|', clusters, fixed = T)]



## filter by KumarMethod
dd.k <- fread(paste0('Cell_Interactions_Sum/Kumar/', data.name, '/interaction_AllBlast_Immune.txt'))
dd.k$clusters = sapply(dd.k$clusters, function(x) gsub('AllBlast', 'Blasts', x) )
pdata.k = NULL
for(cluster0 in unique(pdata$clusters)){
  pdata0 = pdata[clusters == cluster0]
  dd.k0 = dd.k[clusters == cluster0]
  pdata0 = pdata0[pair %in% dd.k0$pair]
  pdata.k = rbind(pdata.k, pdata0)
}


pdata.k$gene1 = sapply(pdata.k$pair, function(x) unlist(strsplit(x, '_'))[1])
pdata.k$gene2 = sapply(pdata.k$pair, function(x) unlist(strsplit(x, '_'))[2])

## filter by DEGs
pdata.k = pdata.k[gene2 %in% degs[cluster == 'Blasts']$gene_name]

## plot overlapped results with kumarMethod
pdata.k$isSig = 1
pdata.mtx = reshape2::recast(pdata.k, pair ~ clusters, measure.var = 'isSig')
rownames(pdata.mtx) = pdata.mtx[, 1]
pdata.mtx = pdata.mtx[, -1]
pdata.mtx[is.na(pdata.mtx)] = 0
pdata.mtx = as.matrix(pdata.mtx)

#pdata.mtx = -log10(pdata.mtx[, -1])
#pdata.mtx[pdata.mtx <= 1] = 1
ph <- pheatmap::pheatmap(pdata.mtx, cluster_cols = T, cluster_rows = T,
                         angle_col = 45, scale = 'none', 
                         clustering_distance_rows = 'euclidean',
                         color = colorRampPalette(c( "white", "#DC143C"))(100))
fname = paste0('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/', data.name, '_relate2_',
               interested.ctype, '_vs_', others.ctype.name, '_heatmap.eps')
ggsave(ph, filename = fname,
       width = 14, height = 7, device = 'eps')

pdata.k[clusters == 'CD8_Cytotoxic|Blasts', 
        'direction' := ifelse(pair %in% names(ann.ptw), ann.ptw[pair], 0)]
pdata.k[clusters != 'CD8_Cytotoxic|Blasts', 
        'direction' := ifelse(pair %in% names(ann.ptw.nk), ann.ptw.nk[pair], 0)]

pdata.k = pdata.k[order(-direction)]
pdata.k$pair = factor(pdata.k$pair, levels = unique(pdata.k$pair))
pdata.k$clusters = factor(pdata.k$clusters, levels = c("CD8_Cytotoxic|Blasts",
                                                       "CD16_NK|Blasts", "CD56_NK|Blasts"))
names(pdata.k)[which(names(pdata.k) == 'mean')] = 'expression'

pd.cyto <- ggplot(pdata.k, aes(x=clusters, y=pair)) +
  geom_point(aes(size=expression, color=expression)) +
  scale_color_gradientn('expression', 
                        colors = viridis(500)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, family = 'Arial'),
        axis.text.y = element_text(size=14, colour = "black", family = 'Arial',
                                   margin = margin(t = 0, r = 0, b = 0, l = 20)),
        axis.title=element_blank(),
        text = element_text('Arial'),
        panel.border = element_rect(size = 0.5, linetype = "solid", colour = "black")) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) 


fname.dot1 = paste0('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/', data.name, '_relate2_',
                    interested.ctype, '_vs_', others.ctype.name,
                    '_dotPlot_cytotoxic.eps')
ggsave(pd.cyto, filename = fname.dot1,
       width = 7, height = 7, device = 'eps')

if(T){
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  myColors = getPalette(15)
  names(myColors) = c('CD14_Monocyte', 'CD16_Monocyte', 'CD16_NK', 'CD4_Early_Activation', 'CD4_Naive', 'CD56_NK', 'CD8_Cytotoxic', 'CD8_EffectorMemory',
                      'CD8_Naive', 'cDC', 'Exhausted_T', 'Mature_B', 'NKT', 'pDC', 'Unknown')
  
  pdata.k$clusters = as.character(pdata.k$clusters)
  pdata.k$ctype1 = sapply(pdata.k$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[1])
  pdata.k$ctype2 = sapply(pdata.k$clusters, function(x) unlist(strsplit(x, '|', fixed = T))[2])
  
  pdata.c <- subset(pdata.k, select = c(gene1, gene2, ctype1, ctype2,
                                                          expression, expression, direction))
  names(pdata.c) = c("receptor", "ligand", "cell_to", "cell_from" , 
                     "cell_to_mean_exprs", "cell_from_mean_exprs", 'direction')
  pdata.c$comm_type = 'checkpoint'
  
  
  cell.col = myColors[unique(pdata.c$cell_to)]
  cell.col['Blasts'] = "red"
  
  pdata.c$link.col = 'gray70'
  pdata.c$link.col = ifelse(pdata.c$direction == -1, 'red',
                            pdata.c$link.col)
  pdata.c$link.col = ifelse(pdata.c$direction == 1, 'skyblue',
                            pdata.c$link.col)
  
  pdata.c = pdata.c[order(direction)]
  postscript('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/circular_plot_BlastvsImmune.eps',
             height = 6, width = 6)
  pdata.c$cell_from = 'Major_Blasts'
  
  cell.col0 = cell.col
  names(cell.col0)[4] = 'Major_Blasts'
  iTALK::LRPlot(pdata.c, datatype = 'mean count', 
                cell_col = cell.col0, transparency = 0,
                link.arr.width = 0.2,
                link.arr.lwd = pdata.c$cell_to_mean_exprs,
                link.arr.col = pdata.c$link.col)
  dev.off()
  
}

## summarize into killing score
if(T){
  pdata.k[, 'kscore' := sum(expression * direction), by = clusters]
  kscore <- pdata.k %>% subset(select = c(clusters, kscore)) %>% .[!duplicated(.)]
  
  kcolor = c('CD8_Cytotoxic|Blasts' = "#C4625D",
             'CD16_NK|Blasts' = "#3A85A8",
             'CD56_NK|Blasts' = "#8D5B96")
  p1 <- ggplot(data = kscore, aes(x = clusters, y = kscore, fill = clusters)) +
    geom_bar(stat = 'identity') + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    NoLegend() + xlab('') + ylab('Pro-killing score') +
    scale_fill_manual(values = kcolor)

  
  write.table(kscore, 
              file = 'Figures/summary/Cell_cell_communication/HSPC1vsBlasts/cytotoxic_pro_killing_score_Blasts.txt',
              row.names = F, quote = F, sep = '\t')
}

## < plot killing score for HSPC1vsBlasts together ####
if(T){
  
  kscore.HSPC1 = fread('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/cytotoxic_pro_killing_score_HSPC1.txt')
  kscore.Blasts = fread('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/cytotoxic_pro_killing_score_Blasts.txt')
  
  kscore.HSPC1$Age = 'HSPC1'
  kscore.Blasts$Age = 'Blasts'
  
  kscore = rbind(kscore.HSPC1, kscore.Blasts)
  kscore$clusters = sapply(kscore$clusters, function(x) gsub('|Blasts', '', x, fixed = T))
  kscore$clusters = sapply(kscore$clusters, function(x) gsub('|HSPC1', '', x, fixed = T))
  
  p2 <- ggplot(data = kscore, aes(x = clusters, y = kscore, fill = Age)) +
    geom_bar(stat = 'identity', position = 'dodge') + theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    xlab('') + ylab('Pro-killing score') +
    scale_fill_manual(values = c("#E41A1C", "#984EA3")) 
  
    #geom_hline(yintercept = 0, linetype = 2)
  
  fname2 = paste0('Figures/summary/Cell_cell_communication/HSPC1vsBlasts/HSPC1vsBlasts_cellphoneDB_Kumar_kscore.eps')
  ggsave(p2, filename = fname2,
         width = 6, height = 6, device = 'eps')
}
