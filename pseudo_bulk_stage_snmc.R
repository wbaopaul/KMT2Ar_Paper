source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)

## list of functions ####
library(entropy)
library(Rcpp)
sourceCpp(code = '
          #include <Rcpp.h>
          using namespace Rcpp;
          // get overlap between two data frame, output a vector, in which 
          // the ith component to be 1 if the ith row of data1 has overlap in data2, zero otherwise
          // [[Rcpp::export]]
          IntegerVector getOverlaps_C1(DataFrame dat1, DataFrame dat2) {
          CharacterVector chr1 = dat1["chr"];
          CharacterVector chr2 = dat2["chr"];
          NumericVector start1 = dat1["start"];
          NumericVector end1 = dat1["end"];
          NumericVector start2 = dat2["start"];
          NumericVector end2 = dat2["end"];
          
          int n1 = chr1.size(), n2 = chr2.size();
          NumericVector midP1(n1), len1(n1), len2(n2), midP2(n2);
          IntegerVector over1(n1);
          
          len1 = (end1 - start1)/2;
          midP1 = (end1 + start1)/2;
          
          len2 = (end2 - start2)/2;
          midP2 = (end2 + start2)/2;
          
          for(int i=0; i<n1; i++){
          over1[i] = 0;
          for(int j=0; j<n2; j++){
          if((chr2[j] == chr1[i]) && (fabs(midP1[i] - midP2[j]) <= max(NumericVector::create(len1[i], len2[j])))){
          over1[i] = 1;
          break;
          }
          }
          }
          
          return(over1);
          }
')


## 0. prepare cell barcodes by stage and age ####
younger_g <-  c('MLLr876533', 'MLLr882304', 'MLLr875706',
                'MLLr1154', 'MLLr870684', 'MLLr879583',
                'MLLr876545', 'MLLr874013', 'MLLr879440',
                'MLLr878289', 'MLLr871427')

final_matching = readRDS('MetaData/snmC/cell_coembed_snmc_gene_body_2kb_referenceATAC.rds')
final_matching$projCtype_merge = final_matching$projCtype
final_matching = final_matching[projCtype %notin% c('NK-like', 'T-like', 'MEP-like')]

final_matching$projCtype_merge[final_matching$projCtype_merge %in% 
                                 c('HSPC-like', 'LMPP-like', 'CLP-like')] = 'Early-Prog-like'
final_matching$projCtype_merge[final_matching$projCtype_merge %in% 
                                 c('cDC-like', 'DC-Progenitor-like', 
                                   'GMP-like', 'Mono-like', 'pDC-like')] = 'M-lineage'
final_matching$lineage = 'B-lineage'
final_matching$lineage[final_matching$projCtype_merge == 'M-lineage'] = 'M-lineage'
final_matching$Age_group = 'Old'
final_matching$Age_group[final_matching$sample %in% younger_g] = 'Young'


## write list of file names into a file
stages = c('Early-Prog-like', 'Pre-pro-B-like',
           'Pro-B-like', 'Pre-B-like', 
           'Immature-B-like', 'Mature-B-like')
groups = c('Young', 'Old')

dir0 = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/snmc/LEUK/pilot/working_dir/Pilot/alignments/'
dir1 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_1/working_dir/'
dir2 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_2/working_dir/'
dir3 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_3/working_dir/'


for(group0 in groups){
  for(stage0 in stages){
    dd = final_matching[Age_group == group0 & projCtype_merge == stage0]
    fpath = sapply(1:nrow(dd), function(x){
      cname = dd[x, snmc_cell]
      sample0 = dd[x, sample]
      ID0 = gsub('MLLr', '', sample0)
      batch0 = dd[x, batch]
      if(batch0 == 'Pilot') return(paste0(dir0, cname, '.organism.sorted.bam'))
      if(batch0 == 'batch1') return(paste0(dir1, ID0, '/alignments/', cname, '.organism.sorted.bam'))
      if(batch0 == 'batch2') return(paste0(dir2, ID0, '/alignments/', cname, '.organism.sorted.bam'))
      if(batch0 == 'batch3') return(paste0(dir3, ID0, '/alignments/', cname, '.organism.sorted.bam'))
    })
    fout = file(paste0('MetaData/snmC/bamfilenames_', group0, '_', stage0, '.txt'))
    writeLines(fpath, fout)
    close(fout)
  }
}

for(group0 in c('B-lineage', 'M-lineage')){
    dd = final_matching[lineage == group0]
    fpath = sapply(1:nrow(dd), function(x){
      cname = dd[x, snmc_cell]
      sample0 = gsub('MLLr', '', dd[x, sample])
      batch0 = dd[x, batch]
      if(batch0 == 'Pilot') return(paste0(dir0, cname, '.organism.sorted.bam'))
      if(batch0 == 'batch1') return(paste0(dir1, sample0, '/alignments/', cname, '.organism.sorted.bam'))
      if(batch0 == 'batch2') return(paste0(dir2, sample0, '/alignments/', cname, '.organism.sorted.bam'))
      if(batch0 == 'batch3') return(paste0(dir3, sample0, '/alignments/', cname, '.organism.sorted.bam'))
    })
    fout = file(paste0('MetaData/snmC/bamfilenames_', group0,  '.txt'))
    writeLines(fpath, fout)
    close(fout)
  
}



## 1. DMR: prepare filtered DMR per stage ####

dir0 = '/mnt/isilon/tan_lab/yuw1/Data/snmc/ByStage/'
stages = c('Early-Prog-like', 'Pre-pro-B-like',
           'Pro-B-like', 'Pre-B-like', 
           'Immature-B-like', 'Mature-B-like')
dmr.young = dmr.old = list()
for(stage0 in stages){
  meth.old = fread(paste0(dir0, 'DMR/DMR_Old_lt_Young_', stage0))
  meth.young = fread(paste0(dir0, 'DMR/DMR_Young_lt_Old_', stage0))
  
  meth.old[, 'V4' := gsub('X:', '', V4), by = V4]
  meth.young[, 'V4' := gsub('X:', '', V4), by = V4]
  names(meth.old)[1:5] = names(meth.young)[1:5] = c('chr', 'start', 'end', 'ntotal', 
                                                    'nless_meth') 
  meth.old$ntotal = as.integer(meth.old$ntotal)
  meth.young$ntotal = as.integer(meth.young$ntotal)
  meth.old[, 'ss' := end -start]
  meth.young[, 'ss' := end -start]
  ## filter
  meth.old = meth.old[ntotal >= 4 & nless_meth >= 2]
  meth.young = meth.young[ntotal >= 4 & nless_meth >= 2]
  dmr.young[[stage0]] = meth.young
  dmr.old[[stage0]] = meth.old
  write.table(meth.young, file = paste0('MetaData/snmC/DMR/dmr_young_', stage0, '.bed'),
              quote = F, row.names = F, sep = '\t')
  write.table(meth.old, file = paste0('MetaData/snmC/DMR/dmr_old_', stage0, '.bed'),
              quote = F, row.names = F, sep = '\t')
  
}



## overlap with DEG by stage ####
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
need.genes = rownames(seurat.rna)
rm(seurat.rna)

gene_ann <- fread("/mnt/isilon/tan_lab/uzuny/projects/sinbad/package/dependencies/annot/hg38/regions.all_genes.bed",
                  select = c(1,2,3, 4, 7, 8))

names(gene_ann) = c('chr', 'start', 'end', 'strand', 'gene_name', 'gene_type')
gene_ann = gene_ann[gene_name %in% need.genes]
gene_ann = gene_ann[gene_type %in% c('protein_coding', 'lincRNA')]
gene_ann = gene_ann[!duplicated(gene_name)]
setkey(gene_ann, gene_name)

## for promoter in gene level
gene_ann[, 'start_prom' := ifelse(strand == '+', start - 2000, end - 2000)]
gene_ann[, 'end_prom' := ifelse(strand == '+', start + 2000, end + 2000)]

## for promoter in transcript level
trx_ann <- fread("/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_tss.bed",
                 select = c(1:4, 6, 7))
names(trx_ann) = c('chr', 'start', 'end',  'gene_name', 'strand', 'gene_type')
trx_ann = trx_ann[gene_type %in% c('protein_coding', 'lincRNA', 'miRNA')]
trx_ann = trx_ann[gene_name %in% need.genes]
trx_ann[, 'Tss' := ifelse(strand == '+', start, end)]
trx_ann[, start := Tss - 2000]
trx_ann[, end := Tss + 2000]

degs.ann.list = list()
deg.dir = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/YoungvsOld/ByProjCtype/'
for(stage0 in stages){
  tmp0 = gsub('-like', '', stage0)
  if(tmp0 == 'Early-Prog') tmp0 = 'EarlyProg'
  degs = data.table(read.table(paste0(deg.dir, '/', tmp0, '_combined_DEG_FC0.5.txt'), header = T))
  names(degs)[1] = 'gene'
  degs$gene = as.character(degs$gene)
  degs$cluster = 'Young'
  degs$cluster[degs$avg_logFC < 0 ] = 'Old'
  #degs = degs[cluster == group0]
  degs = degs[order(-abs(avg_logFC))]
  
  degs = degs[abs(avg_logFC) > 0.5 & p_val_adj < 0.01]
  
  ## add gene annotation
  degs[, 'chr' := gene_ann[J(degs$gene)]$chr]
  degs[, 'start_prom' := gene_ann[J(degs$gene)]$start_prom]
  degs[, 'end_prom' := gene_ann[J(degs$gene)]$end_prom]
  
  degs[, 'start' := gene_ann[J(degs$gene)]$start]
  degs[, 'end' := gene_ann[J(degs$gene)]$end]
  
  meth.old = dmr.old[[stage0]]
  meth.young = dmr.young[[stage0]]
  
  degs$genebody.overlap.old.dmethy = getOverlaps_C1(degs, meth.old)
  degs$genebody.overlap.young.dmethy = getOverlaps_C1(degs, meth.young)
  
  ## overlap with promoter -- gene level
  degs0 = subset(degs, select = c(chr, start_prom, end_prom, gene))
  names(degs0)[2:3] = c('start', 'end')
  degs$promoter.overlap.old.dmethy = getOverlaps_C1(degs0, meth.old)
  degs$promoter.overlap.young.dmethy = getOverlaps_C1(degs0, meth.young)
  
  ## overlap with promoter -- transcript level
  trx.young.dmethy = trx.old.dmethy = rep(0, nrow(degs0))
  for(i in 1:nrow(degs0)){
    chr0 = degs0[i]$chr
    gene0 = degs0[i]$gene
    trx0 = trx_ann[gene_name == gene0]
    if(nrow(trx0) == 0) next
    meth.young0 = meth.young[chr == chr0]
    meth.old0 = meth.old[chr == chr0]
    over.young = getOverlaps_C1(trx0, meth.young0)
    over.old = getOverlaps_C1(trx0, meth.old0)
    if(any(over.young == 1)) trx.young.dmethy[i] = 1
    if(any(over.old == 1)) trx.old.dmethy[i] = 1
  }
  degs$trx.promoter.overlap.young.dmethy = trx.young.dmethy
  degs$trx.promoter.overlap.old.dmethy = trx.old.dmethy

  degs$stage = stage0
  degs.ann.list[[stage0]] = degs
}

degs.comb = do.call('rbind', degs.ann.list)

degs.comb.young = degs.comb[cluster == 'Young']
degs.comb.young.dmr = degs.comb.young[promoter.overlap.young.dmethy == 1 | genebody.overlap.young.dmethy == 1]

degs.comb.old = degs.comb[cluster == 'Old']
degs.comb.old.dmr = degs.comb.old[promoter.overlap.old.dmethy == 1 | genebody.overlap.old.dmethy == 1]

saveRDS(degs.comb, file = 'MetaData/snmC/Summary/degs_with_dmr.rds')
saveRDS(degs.comb.young, file = 'MetaData/snmC/Summary/degs_upYoung_with_dmr.rds')
saveRDS(degs.comb.old, file = 'MetaData/snmC/Summary/degs_upOld_with_dmr.rds')

## < plot for promoter ####
fracs.young.promoter <- sapply(stages, function(x) {
  nn <- nrow(degs.comb.young[trx.promoter.overlap.young.dmethy == 1 & stage == x])
  nd <- nrow(degs.comb.young[ stage == x])
  return(nn/nd)
})

fracs.old.promoter <- sapply(stages, function(x) {
  nn <- nrow(degs.comb.old[trx.promoter.overlap.old.dmethy == 1 & stage == x])
  nd <- nrow(degs.comb.old[ stage == x])
  return(nn/nd)
})

data.promoter = data.table('stage' = rep(names(fracs.young.promoter), 2),
                           'frac' = c(fracs.young.promoter, fracs.old.promoter),
                           'direction' = rep(c('up_young', 'down_young'), each = 6))

data.promoter$stage = factor(data.promoter$stage, levels = names(fracs.old.promoter))
p1 <- ggplot(data = data.promoter[direction == 'up_young'], 
       aes(x = stage, y = frac, fill = direction)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  ggtitle('% Promoter less methylated for up-DEGs in Young') + 
  scale_fill_manual(values = "#E41A1C") +
  geom_text(aes(label=paste0(round(frac*100, 2), '%')), vjust=-0,  
            size=3.5, hjust = 0.5) +
  theme_classic() + NoLegend()

p2 <- ggplot(data = data.promoter[direction == 'down_young'], 
       aes(x = stage, y = frac, fill = direction)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  ggtitle('% Promoter more methylated for down-DEGs in Young') +
  scale_fill_manual(values = "#377EB8") +
  geom_text(aes(label=paste0(round(frac*100, 2), '%')), vjust=-0,  
            size=3.5, hjust = 0.5) +
  theme_classic() + NoLegend()

pcomb <- ggpubr::ggarrange(p1, p2, ncol = 1)
ggsave(pcomb, filename = 'Figures/summary/snmC/frac_degs_promoter_demethylated.eps',
       width = 7, height = 5)

degs.promoter.dmeth <- rbind(degs.comb.young[trx.promoter.overlap.young.dmethy == 1],
                             degs.comb.old[trx.promoter.overlap.old.dmethy == 1])

mtx.promoter = subset(degs.promoter.dmeth,
                      select = c(stage, gene, cluster))

mtx.promoter[, 'direction' := ifelse(cluster == 'Young', 1, -1)]

df.promoter <- data.table::dcast(data = mtx.promoter, stage~gene, fill = 0,
                                 value.var = 'direction')

rnames = df.promoter$stage

df.promoter1 = as.matrix(df.promoter[, -1])

rownames(df.promoter1) = rnames
df.promoter1 = df.promoter1[c('Early-Prog-like', 'Pre-pro-B-like',
                              'Pro-B-like', 'Pre-B-like', 'Immature-B-like'), ]

ann_column = data.frame('stage' = stages,
                        'barcode' = stages,
                        stringsAsFactors = F)
rownames(ann_column) = stages
ann_column$barcode <- NULL
ann_color = list('stage' = c("Early-Prog-like" = "#E3191C", "Pre-pro-B-like" = "#347EB5",
                             "Pro-B-like" = "#4CAC47", "Pre-B-like" = "#A35728", 
                             "Immature-B-like" = "#924C9F", "Mature-B-like" = "#00FFFF"))


ph1 <- pheatmap::pheatmap(as.matrix(t(df.promoter1)), cluster_cols = F,
                   border_color = "NA", annotation_col = ann_column,
                   annotation_colors = ann_color, show_colnames = F,
                   colorRampPalette(c("blue", "white", "red"))(3))
ggsave(ph1, filename = 'Figures/summary/snmC/heatmap_degs_promoter_demethylated.eps',
       width = 8, height = 10)


## < plot for gene body ####

fracs.young.genebody <- sapply(stages, function(x) {
  nn <- nrow(degs.comb.young[genebody.overlap.young.dmethy == 1 & stage == x])
  nd <- nrow(degs.comb.young[ stage == x])
  return(nn/nd)
})

fracs.old.genebody <- sapply(stages, function(x) {
  nn <- nrow(degs.comb.old[genebody.overlap.old.dmethy == 1 & stage == x])
  nd <- nrow(degs.comb.old[ stage == x])
  return(nn/nd)
})


data.genebody = data.table('stage' = rep(names(fracs.young.genebody), 2),
                           'frac' = c(fracs.young.genebody, fracs.old.genebody),
                           'direction' = rep(c('up_young', 'down_young'), each = 6))

data.genebody$stage = factor(data.genebody$stage, levels = names(fracs.old.genebody))
p1 <- ggplot(data = data.genebody[direction == 'up_young'], 
             aes(x = stage, y = frac, fill = direction)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  ggtitle('% genebody more methylated for up-DEGs in Young') + 
  scale_fill_manual(values = "#E41A1C") +
  geom_text(aes(label=paste0(round(frac*100, 2), '%')), vjust=-0,  
            size=3.5, hjust = 0.5) +
  theme_classic() + NoLegend()

p2 <- ggplot(data = data.genebody[direction == 'down_young'], 
             aes(x = stage, y = frac, fill = direction)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  ggtitle('% genebody less methylated for down-DEGs in Young') +
  scale_fill_manual(values = "#377EB8") +
  geom_text(aes(label=paste0(round(frac*100, 2), '%')), vjust=-0,  
            size=3.5, hjust = 0.5) +
  theme_classic() + NoLegend()

pcomb <- ggpubr::ggarrange(p1, p2, ncol = 1)
ggsave(pcomb, filename = 'Figures/summary/snmC/frac_degs_genebody_demethylated.eps',
       width = 7, height = 5)

degs.genebody.dmeth <- rbind(degs.comb.young[genebody.overlap.young.dmethy == 1],
                             degs.comb.old[genebody.overlap.old.dmethy == 1])

mtx.genebody = subset(degs.genebody.dmeth,
                      select = c(stage, gene, cluster))

mtx.genebody[, 'direction' := ifelse(cluster == 'Young', 1, -1)]

df.genebody <- data.table::dcast(data = mtx.genebody, stage~gene, fill = 0,
                                 value.var = 'direction')

rnames = df.genebody$stage

df.genebody1 = as.matrix(df.genebody[, -1])

rownames(df.genebody1) = rnames
df.genebody1 = df.genebody1[c('Early-Prog-like', 'Pre-pro-B-like',
                              'Pro-B-like', 'Pre-B-like', 'Immature-B-like'), ]

ann_column = data.frame('stage' = stages,
                        'barcode' = stages,
                        stringsAsFactors = F)
rownames(ann_column) = stages
ann_column$barcode <- NULL
ann_color = list('stage' = c("Early-Prog-like" = "#E3191C", "Pre-pro-B-like" = "#347EB5",
                             "Pro-B-like" = "#4CAC47", "Pre-B-like" = "#A35728", 
                             "Immature-B-like" = "#924C9F", "Mature-B-like" = "#00FFFF"))


ph2 <- pheatmap::pheatmap(as.matrix(t(df.genebody1)), cluster_cols = F,
                   border_color = "NA", annotation_col = ann_column,
                   annotation_colors = ann_color, show_colnames = F,
                   colorRampPalette(c("blue", "white", "red"))(3))
ggsave(ph2, filename = 'Figures/summary/snmC/heatmap_degs_genebody_demethylated.eps',
       width = 8, height = 12)





## DMR overlap with enhancer by stage specific ep -- not used ####
motif_pk_match = readRDS('EP_Prediction/motif_pk_match_mtx.rds')

ep.young.list = ep.old.list = list()
for(stage0 in stages){
  ep.young = fread('EP_Prediction/regrRes4_EP_in_Young_group.txt')
  ep.old = fread('EP_Prediction/regrRes4_EP_in_Old_group.txt')
  ## sele peaks
  
  final.peaks.young = readRDS(paste0('MetaData/scATAC/peaks_Young_', 
                                     stage0, '_meanUpFC1.414_YoungvsOld.rds'))
  final.peaks.young = sapply(final.peaks.young, function(x) unlist(strsplit(x, ','))[1])
  names(final.peaks.young) = NULL
  ep.young = ep.young[enhancer_peak %in% final.peaks.young]
  
  final.peaks.old = readRDS(paste0('MetaData/scATAC/peaks_Old_',
                                   stage0, '_meanUpFC1.414_YoungvsOld.rds'))
  final.peaks.old = sapply(final.peaks.old, function(x) unlist(strsplit(x, ','))[1])
  names(final.peaks.old) = NULL
  ep.old = ep.old[enhancer_peak %in% final.peaks.old]
  
  ## select degs
  deg.dir = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/YoungvsOld/ByProjCtype/'
  tmp0 = gsub('-like', '', stage0)
  if(tmp0 == 'Early-Prog') tmp0 = 'EarlyProg'
  degs = data.table(read.table(paste0(deg.dir, '/', tmp0, '_combined_DEG_FC0.5.txt'), header = T))
  names(degs)[1] = 'gene'
  degs$cluster = 'Young'
  degs$cluster[degs$avg_logFC < 0 ] = 'Old'
  
  degs = degs[order(-abs(avg_logFC))]
  
  ep.young = ep.young[gene_name %in% degs[cluster == 'Young']$gene]
  ep.old = ep.old[gene_name %in% degs[cluster == 'Old']$gene]
  
  
  ## < filter by TF hits at enhancer side 
  
  enriched.tf.list = readRDS('EP_Prediction/Blasts_YoungvsOld_enriched_tfs_list_by_stage_exprfrac0.1.rds')
  enriched.tf1 = enriched.tf.list[[stage0]]
  enriched.tf1.young = enriched.tf1[cluster1 == 'Young', ]
  enriched.tf1.old = enriched.tf1[cluster1 == 'Old', ]
  
  meth.young = dmr.young[[stage0]]
  meth.old = dmr.old[[stage0]]
  
  if(nrow(enriched.tf1.young) != 0) {
    sele.tf.mat = motif_pk_match[, enriched.tf1.young$feature] 
    sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
    ep.young = ep.young[enhancer_peak %in% sele.peaks]
    ep.young = subset(ep.young, select = c(gene_name, enhancer_peak, ep_dist, fdr, Estimate))
    ep.young$stage = stage0
    
    ep.loc = tidyr::separate(data = subset(ep.young, select = enhancer_peak), 
                             into = c('chr', 'start', 'end'), col = 'enhancer_peak')
    ep.loc$start = as.integer(ep.loc$start)
    ep.loc$end = as.integer(ep.loc$end)
    ep.young$young.dmethy.overlap.enh = getOverlaps_C1(ep.loc, meth.young)
    ep.young$old.demth.overlap.enh = getOverlaps_C1(ep.loc, meth.old)
    
    ep.young.list[[stage0]] = ep.young
    
  }
  
  if(nrow(enriched.tf1.old) != 0) {
    sele.tf.mat = motif_pk_match[, enriched.tf1.old$feature] 
    sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
    ep.old = subset(ep.old, select = c(gene_name, enhancer_peak, ep_dist, fdr, Estimate))
    ep.old = ep.old[enhancer_peak %in% sele.peaks]
    ep.old$stage = stage0
    ep.loc = tidyr::separate(data = subset(ep.old, select = enhancer_peak), 
                             into = c('chr', 'start', 'end'), col = 'enhancer_peak')
    ep.loc$start = as.integer(ep.loc$start)
    ep.loc$end = as.integer(ep.loc$end)
    ep.old$old.dmethy.overlap.enh = getOverlaps_C1(ep.loc, meth.old)
    ep.old$young.dmethy.overlap.enh = getOverlaps_C1(ep.loc, meth.young)
    ep.old.list[[stage0]] = ep.old
  }
  
  
  
  
}

ep.young.comb = do.call('rbind', ep.young.list)
ep.old.comb = do.call('rbind', ep.old.list)

saveRDS(ep.young.comb, file = 'MetaData/snmC/Summary/ep_upYoung_with_dmr_atEnhancer.rds')
saveRDS(ep.old.comb, file = 'MetaData/snmC/Summary/ep_upOld_with_dmr_atEnhancer.rds')







## DMR overlap with enhancer enriched with stage specific TF ####
motif_pk_match = readRDS('EP_Prediction/motif_pk_match_mtx.rds')

enh.young.list = enh.old.list = list()
for(stage0 in stages){
  ## sele peaks
  
  final.peaks.young = readRDS(paste0('MetaData/scATAC/peaks_Young_', 
                                     stage0, '_meanUpFC1.414_YoungvsOld.rds'))
  final.peaks.young = sapply(final.peaks.young, function(x) unlist(strsplit(x, ','))[1])
  names(final.peaks.young) = NULL
  
  final.peaks.old = readRDS(paste0('MetaData/scATAC/peaks_Old_',
                                   stage0, '_meanUpFC1.414_YoungvsOld.rds'))
  final.peaks.old = sapply(final.peaks.old, function(x) unlist(strsplit(x, ','))[1])
  names(final.peaks.old) = NULL
  
  ## < filter by TF hits at enhancer side 
  
  enriched.tf.list = readRDS('EP_Prediction/Blasts_YoungvsOld_enriched_tfs_list_by_stage_exprfrac0.1.rds')
  enriched.tf1 = enriched.tf.list[[stage0]]
  enriched.tf1.young = enriched.tf1[cluster1 == 'Young', ]
  enriched.tf1.old = enriched.tf1[cluster1 == 'Old', ]
  
  meth.young = dmr.young[[stage0]]
  meth.old = dmr.old[[stage0]]
  
  if(nrow(enriched.tf1.young) != 0) {
    sele.tf.mat = motif_pk_match[, enriched.tf1.young$feature] 
    sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
    
    sele.tf.mat = sele.tf.mat[sele.peaks, ]
    
    peak.young = final.peaks.young[final.peaks.young %in% sele.peaks]
    sele.tf.mat = sele.tf.mat[peak.young, ]
    
    enh.young = tidyr::separate(data = data.table('enhancer_peak' = peak.young), 
                                 into = c('chr', 'start', 'end'), 
                                col = 'enhancer_peak', remove = F)
    enh.young$start = as.integer(enh.young$start)
    enh.young$end = as.integer(enh.young$end)
    enh.young$stage = stage0
    
    enh.young$overlap.young.dmethy = getOverlaps_C1(enh.young, meth.young)
    enh.young$overlap.old.dmethy = getOverlaps_C1(enh.young, meth.old)
    
    enh.young.list[[stage0]] = enh.young
    
  }
  
  if(nrow(enriched.tf1.old) != 0) {
    sele.tf.mat = motif_pk_match[, enriched.tf1.old$feature] 
    sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
    peak.old = final.peaks.old[final.peaks.old %in% sele.peaks]
    enh.old = tidyr::separate(data = data.table('enhancer_peak' = peak.old), 
                                into = c('chr', 'start', 'end'), 
                              col = 'enhancer_peak', remove = F)
    enh.old$start = as.integer(enh.old$start)
    enh.old$end = as.integer(enh.old$end)
    enh.old$stage = stage0
    
    enh.old$overlap.old.dmethy = getOverlaps_C1(enh.old, meth.old)
    enh.old$overlap.old.dmethy = getOverlaps_C1(enh.old, meth.young)
    enh.old.list[[stage0]] = enh.old
    
  }
  
  
  
  
}

enh.young.comb = do.call('rbind', enh.young.list)
enh.old.comb = do.call('rbind', enh.old.list)

enh.young.comb = enh.young.comb[overlap.young.dmethy == 1]
enh.old.comb = enh.old.comb[overlap.old.dmethy == 1]

saveRDS(enh.young.comb, file = 'MetaData/snmC/Summary/enh_upYoung_with_dmr.rds')
saveRDS(enh.old.comb, file = 'MetaData/snmC/Summary/enh_upOld_with_dmr.rds')

lapply(enh.young.list, nrow)
lapply(enh.old.list, nrow)

lapply(enh.young.list, function(x) sum(x$overlap.young.dmethy))
lapply(enh.old.list, function(x) sum(x$overlap.old.dmethy))

enriched.tf.list = readRDS('EP_Prediction/Blasts_YoungvsOld_enriched_tfs_list_by_stage_exprfrac0.1.rds')
enriched.tf.comb = lapply(names(enriched.tf.list), function(x) {
  enriched.tf.list[[x]]$stage = x
  return(enriched.tf.list[[x]])
  }
)
enriched.tf.comb = do.call('rbind', enriched.tf.comb)

## construct stage-by-TF table, record the number of enhancers with TF hits and
## overlapped with demethlyated region in young/old group, resp
stage.tf.dmethy.young = stage.tf.dmethy.old = list()
for(stage0 in stages){
  dd.young = enh.young.comb[stage == stage0]
  dd.old = enh.old.comb[stage == stage0]
  
  tfs.young = enriched.tf.comb[stage == stage0 & cluster1 == 'Young']$feature
  tfs.old = enriched.tf.comb[stage == stage0 & cluster1 == 'Old']$feature
  
  tf.peak.young = motif_pk_match[dd.young$enhancer_peak, tfs.young]
  tf.peak.old = motif_pk_match[dd.old$enhancer_peak, tfs.old]
  
  tf.dmethy.young = colSums(tf.peak.young)
  tf.dmethy.old = colSums(tf.peak.old)
  
  if(length(tf.dmethy.young) > 0) tf.dmethy.young = data.table('stage' = stage0,
                                                             'TF' = names(tf.dmethy.young),
                                                             'N' = as.vector(tf.dmethy.young))
  
  if(length(tf.dmethy.old) > 0) tf.dmethy.old = data.table('stage' = stage0,
                                                             'TF' = names(tf.dmethy.old),
                                                             'N' = -as.vector(tf.dmethy.old))
  
  
  stage.tf.dmethy.young[[stage0]] = tf.dmethy.young
  stage.tf.dmethy.old[[stage0]] = tf.dmethy.old
  
  
}
stage.tf.dmethy.young[['Early-Prog-like']] = NULL
stage.tf.dmethy = rbind(do.call('rbind', stage.tf.dmethy.young),
                        do.call('rbind', stage.tf.dmethy.old))

stage.dmethy = dcast(data = stage.tf.dmethy, stage ~ TF, value.var = 'N', fill = 0)
rnames = stage.dmethy$stage
stage.dmethy = as.matrix(stage.dmethy[, -1])
rownames(stage.dmethy) = rnames
stage.dmethy = stage.dmethy[stages, ]

ann_column = data.frame('stage' = stages,
                        'barcode' = stages,
                        stringsAsFactors = F)
rownames(ann_column) = stages
ann_column$barcode <- NULL
ann_color = list('stage' = c("Early-Prog-like" = "#E3191C", "Pre-pro-B-like" = "#347EB5",
                             "Pro-B-like" = "#4CAC47", "Pre-B-like" = "#A35728", 
                             "Immature-B-like" = "#924C9F", "Mature-B-like" = "#00FFFF"))


ann_color$stage = ann_color$stage[-6]
ph2 <- pheatmap::pheatmap(t(stage.dmethy[-6, ]), cluster_cols = F,
                          border_color = "NA", annotation_col = ann_column[-6],
                          annotation_colors = ann_color, show_colnames = F,
                          breaks = sort(unique(c(seq(-400, 0, length.out = 50), 
                                                 seq(0, 400, length.out = 50)))),
                          colorRampPalette(c("blue", "white", "orange"))(99))
ggsave(ph2, filename = 'Figures/summary/snmC/heatmap_TF_enhancerNumb_overlap_demethylated.eps',
       width = 8, height = 10)


