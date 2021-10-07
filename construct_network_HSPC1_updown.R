source('scDataAnalysis_Utilities.R')
library(igraph)
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

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



## construct TRN for a specify condition ####
motif_pk_match = readRDS('EP_Prediction/motif_pk_match_mtx.rds')


loop_type = 'loops_filtered_from_overall'  

  cell_type = 'HPSC1'
  predicted.ep = fread('EP_Prediction/regrRes4_EP_overall.txt')
  
  ## sele peaks
  
  final.peaks.str = readRDS('MetaData/scATAC/peaks_hspc1_meanUpFC1.414_HSPCvsHSPC1.rds')
  dt = tidyr::separate(predicted.ep, col = 'enhancer_peak', remove = F,
                       into = c("chr", "start", "end")) %>% subset(select = c(chr, start, end))
  
  final.peaks.str = tidyr::separate(data.table(x = final.peaks.str), col = 'x',
                                    into = c("chr", "start", "end"))
  class(dt$start) = class(dt$end) = 'integer'
  class(final.peaks.str$start) = class(final.peaks.str$end) = 'integer'
  overlap.res = getOverlaps_C1(dt, final.peaks.str)
  predicted.ep = predicted.ep[overlap.res == 1]
  
  
  ## select degs
 
  degs = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/stagewise_DEG_5HDProjection/0_HSPC1_DEG_HSPC_Progenitors_LR_Alloutput.txt')
  names(degs)[1] = 'gene'
  degs = degs[, 'cluster' := ifelse(avg_logFC > 0, 'HSPC', 'HSPC1')]
  degs1 = degs[avg_logFC < -0.5 & p_val_adj < 0.05]
  degs1 = degs1[order(p_val_adj)]
  degs1 = degs1[1:30, ]
  degs2 = degs[avg_logFC > 0.5 & p_val_adj < 0.05]
  
  degs = rbind(degs1, degs2)
  
  predicted.ep = predicted.ep[gene_name %in% degs$gene]
  
  
  ## < filter by TF hits at enhancer side ####
  
  enriched.tf1 = readRDS('EP_Prediction/HSPC1_enriched_tfs_HSPCvsHSPC1.rds')[[1]]
  #enriched.tf = enriched.tf1$feature[1:30]
  
  sele.tfs = c('JUN', 'JUND', 'JUNB', 'FOS', 'FOSB',
               "IRF1", "IRF3", "IRF7", "IRF9",
               "STAT1", "STAT3", 'STAT5A', 'STAT5B',
               'NFKB1', 'NFKB2', 'REL', 'RELA', 'RELB')
  enriched.tf1 = enriched.tf1[feature %in% sele.tfs, ]
  enriched.tf = enriched.tf1$feature
  
  enriched.tf = enriched.tf[!is.na(enriched.tf)]
  sele.tf.mat = motif_pk_match[, enriched.tf] 
  sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
  
  predicted.ep = predicted.ep[enhancer_peak %in% sele.peaks]
  
  write.table(predicted.ep, file = paste0('EP_Prediction/EP4UCSC/', cell_type, '_regr_ep_updown_seleTFs.txt'),
              row.names = F, quote = F, sep = '\t')
  
  
  if(nrow(predicted.ep) == 0) next
  ## < construct network ####
  ## split loop by tf
  ep.tf = list()
  for(TF0 in enriched.tf){
    peaks0 = names(which(sele.tf.mat[, TF0] > 0))
    ep0 = predicted.ep[enhancer_peak %in% peaks0]
    ep0 = subset(ep0, select = c('gene_name', 'Estimate', 'fdr'))
    ep0[, 'score' := -log10(fdr)]
    ep0$TF = TF0
    ep0[, 'N' := .N, by = gene_name]
    ep0[, 'score' := sum(score), by = gene_name]
    ep0 = subset(ep0, select = c('TF', 'gene_name', 'score'))
    ep.tf[[TF0]] = ep0[!duplicated(ep0)]
  }
  
  ep.tf.comb = do.call('rbind', ep.tf)
  
  ep.tf.comb = ep.tf.comb[order(-score)]
  ep.tf.comb = ep.tf.comb[1:100, ]
  
  vertex.gr = data.table('gene_name' = unique(c(ep.tf.comb$TF, ep.tf.comb$gene_name)),
                         'group' = 'Gene')
  vertex.gr$group[vertex.gr$gene_name %in% enriched.tf] = 'TF'
  setkey(vertex.gr, gene_name)
  
  ## add expression information
  
  load('MetaData/scRNA/efreq_ctypte_mllr.RData')
  efreq0 <- efreq.hspc1
  
  
  enriched.tf1[, 'delta' := mean1-mean0]
  
  vertex.gr[gene_name %in% degs$gene, 
            'expr_logFC' := -(degs[gene == gene_name]$avg_logFC), by = gene_name]
  vertex.gr[gene_name %in% degs$gene, 
            'log10PV' := -log10(degs[gene == gene_name]$p_val_adj), by = gene_name]
  
  vertex.gr[gene_name %in% enriched.tf, 
            'tf_dev' := abs(enriched.tf1[feature == gene_name]$delta), by = gene_name]
  vertex.gr[gene_name %in% enriched.tf, 
            'log10PV' := -log10(enriched.tf1[feature == gene_name]$pv_adjust), by = gene_name]
  
  vertex.gr$effect_size = vertex.gr$tf_dev /max(vertex.gr$tf_dev, na.rm = T)
  vertex.gr$tmp = vertex.gr$expr_logFC /max(abs(vertex.gr$expr_logFC), na.rm = T)
  vertex.gr[, 'effect_size' := ifelse(is.na(tf_dev), tmp, effect_size), by = gene_name]
  vertex.gr[, 'tmp' := NULL]
  vertex.gr[is.infinite(log10PV)]$log10PV = max(vertex.gr[!is.infinite(log10PV)]$log10PV)+10
  
  ep.tf.comb$direction = 'up'
  ep.tf.comb[gene_name %in% degs2$gene]$direction = 'down'
  
  #ep.tf.comb[, 'score' := ifelse(direction == 'down', -score, score)]
  
  vertex.gr$direction = 'up'
  vertex.gr[gene_name %in% degs2$gene & group == 'Gene']$direction = 'down'
  
  vertex.gr[, 'score' := ifelse(direction == 'down', -log10PV, log10PV)]
  
  write.table(ep.tf.comb, file = paste0('EP_Prediction/GRN/', loop_type, '/edges4grn_', cell_type, '_updown_seleTFs_top100EP.txt'),
              sep = '\t', row.names = F, quote = F)
  write.table(vertex.gr, file = paste0('EP_Prediction/GRN/', loop_type, '/vertices4grn_', cell_type, '_updown_seleTFs_top100EP.txt'),
              sep = '\t', row.names = F, quote = F)
  

  
  ## write as a supplementary table (using all in the list)
  setkey(enriched.tf1, feature)
  final.ep.tf.table = ep.tf.comb
  final.ep.tf.table[, 'TF_diff_deviation' := enriched.tf1[feature == TF]$delta, by = TF]
  final.ep.tf.table[, 'TF_FDR' := enriched.tf1[feature == TF]$pv_adjust, by = TF]
  final.ep.tf.table[, 'gene_log2FC' := degs[gene == gene_name]$avg_logFC, 
                    by = gene_name]
  final.ep.tf.table[, 'gene_FDR' := degs[gene == gene_name]$p_val_adj, 
                    by = gene_name]
  names(final.ep.tf.table)[3] = '-log10FDR(EP)'
  
  write.table(final.ep.tf.table, file = paste0('EP_Prediction/GRN/', loop_type, 
                                               '/TRNtable_HSPC1vsHSPC.txt'),
              sep = '\t', row.names = F, quote = F)
