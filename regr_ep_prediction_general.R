source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)
## get a binary matrix indicates the gene-peak affinity
## gene.list are data.table including column name gene_name 
## gene_ann should include gene_name,chr,start,end
get_gene2peak_map <- function(gene.list, peak_names,
                              gene_ann, distal_dist = 2e05){
  
  peaks <- tidyr::separate(data.table(peak_name = peak_names), 
                           col = peak_name, remove = F,
                           into = c('chr', 'start', 'end'))
  class(peaks$start) = class(peaks$end) = 'integer'
  setkey(gene_ann, gene_name)
  setkey(peaks, peak_name)
  gene.list = gene.list[gene_name %in% gene_ann$gene_name]
  gene.list$chr = gene_ann[gene.list$gene_name]$chr
  gene.list$start = gene_ann[gene.list$gene_name]$start
  gene.list$end = gene_ann[gene.list$gene_name]$end
  
  ## for each gene, get the corresponding peaks
  gene2peaks = lapply(gene.list$gene_name, function(x) {
    
    chr0 = gene_ann[x]$chr
    start0 = gene_ann[x]$start
    end0 = gene_ann[x]$end
    
    peaks0 = peaks[chr == chr0]
    peaks0 = peaks0[abs(start/2 + end/2 - start0/2 - end0/2) <= distal_dist]
    return(peaks0$peak_name)
  } )
  
  ## pool all peaks relate to one gene
  gene2peaks.u <- lapply(sort(unique(gene.list$gene_name)), function(x){
    id = which(gene.list$gene_name == x)
    tmp_peak <- do.call('c', lapply(id, function(x) gene2peaks[[x]]))
    return(tmp_peak)
  })
  names(gene2peaks.u) <- sort(unique(gene.list$gene_name))
  lens = sapply(gene2peaks.u, length)
  
  genes.f <- names(which(lens > 0))
  lens = lens[lens > 0]
  ## construct overlap matrix
  gene2peaks.dt <- data.table('gene' = rep(genes.f, lens),
                              'peak' = do.call('c', lapply(genes.f, 
                                                           function(x) gene2peaks.u[[x]])))
  upeaks = sort(unique(gene2peaks.dt$peak))
  gene2peaks.dt[, 'id1' := which(genes.f == gene), by = gene]
  gene2peaks.dt[, 'id2' := which(upeaks == peak), by = peak]
  gene2peak.map <- sparseMatrix(i = gene2peaks.dt$id1,
                                j = gene2peaks.dt$id2,
                                dimnames = list(genes.f, upeaks))
  gene2peak.map = gene2peak.map * 1
  
  return(gene2peak.map)
}

## annotate peaks with gene +/- 5kb of its TSS
# input peak_coords with chr-start-end, format
annPeak2Gene <- function(peak_coords, gene_ann, proxim_dist = 5e+03){
  gene_ann[, 'tss' := ifelse(strand == '+', start, end)]
  peaks = tidyr::separate(data.table(x = peak_coords),
                          col = x,
                          into = c('chr', 'start', 'end'))
  
  
  peaks$peak_name = peak_coords
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(peaks$chr)
  peaks_ann = NULL
  for(chr0 in chrs){
    peaks0 = peaks[chr == chr0]
    genes0 = gene_ann[chr == chr0]
    
    peaks0$gene_name = ''
    for(i in 1:nrow(peaks0)){
      tss0 = genes0[tss <= (peaks0$end[i] + proxim_dist) & tss >= (peaks0$start[i] - proxim_dist)]
      if(nrow(tss0) > 0 ) {
        peaks0$gene_name[i] = paste(unique(tss0$gene_name), collapse = ',')
      }
    }
    
    peaks_ann = rbind(peaks_ann, peaks0)
  }
  
  peaks_ann[, 'peak_new_name' := ifelse(!is.na(gene_name) & nchar(gene_name) > 1, 
                                        paste0(peak_name, ',', gene_name), peak_name)]
  
  
  setkey(peaks_ann, peak_name)
  
  return(peaks_ann)
  
}


## map gene to overlapping atac peak
## gene_list with genename, chr, start, end
geneOverlapPeak <- function(gene_list, peak_names, mid_dist = 1000){
  # should include tss information in gene_list
  peaks = tidyr::separate(data = data.table('peak_name' = peak_names),
                          col = peak_name, into = c('chr', 'start', 'end'),
                          remove = F)
  class(peaks$chr) = 'character'
  class(peaks$start) = 'integer'
  class(peaks$end) = 'integer'
  
  chrs = unique(gene_list$chr)
  gene_new = NULL
  peaks[, 'midP' := start/2 + end/2]
  for(chr0 in chrs){
    gene0 = gene_list[chr == chr0, ]
    gene0$peak_name = 'Not_Found'
    peaks0 = peaks[chr == chr0]
    gene0[, 'peak_id0' := any( abs(peaks0$midP -start) < mid_dist | abs(peaks0$midP - end) < mid_dist), 
          by = gene_name]
    gene1 = gene0[peak_id0 == FALSE]
    gene2 = gene0[peak_id0 == TRUE]
    gene2[, 'peak_id' :=  which.min(abs(peaks0$midP - start - 1000)), by = gene_name]
    
    gene2[, 'peak_name' :=  peaks0[peak_id]$peak_name, by = gene_name]
    gene2$peak_id = NULL
    gene_new = rbind(gene_new, gene1, gene2)
  }
  gene_new[, c('peak_id0') := NULL]
  return(gene_new)
}


## seurat co-embedding ####
seurat.rna <- readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')

seurat.atac <- readRDS('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap10000.rds')

## downsample 40K cells --rna
set.seed(2020)
seurat.rna$bc = colnames(seurat.rna)
sele.cells = sample(seurat.rna$bc, 40000)
seurat.rna = subset(seurat.rna, bc %in% sele.cells)
seurat.rna$bc <- NULL


## downsample 40K cells --atac
set.seed(2019)
seurat.atac$bc = colnames(seurat.atac)
sele.cells = sample(seurat.atac$bc, 40000)
seurat.atac = subset(seurat.atac, bc %in% sele.cells)
seurat.atac$bc <- NULL


atac.mtx = seurat.atac@assays$ATAC@counts
rn = rownames(atac.mtx)
rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/yuw1/local_tools/annotation/GRCh38_genes.gtf',
                                            atac.mtx, 
                                            include_body = T)

seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)
genes4anchors = VariableFeatures(object = seurat.rna)



DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- FindVariableFeatures(seurat.atac)
seurat.atac <- ScaleData(seurat.atac)

DefaultAssay(seurat.atac) <- "ATAC"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

## transfer label 

#genes4anchors = NULL
transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)

#co-embedding
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.rna, assay = "RNA", slot = "data")[genes4anchors, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.atac[["pca"]],
                           dims = 1:ncol(seurat.atac[["pca"]]))

# this line adds the imputed data matrix to the seurat.atac object
seurat.atac[["RNA"]] <- imputation
coembed <- merge(x = seurat.rna, y = seurat.atac)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
DimPlot(coembed, group.by = 'projCtype', label = T) 
saveRDS(coembed, file = 'Seurat_Objects/Integrated/seurat_18MLLr_40Kcoembed.rds')

## 1 to 1 cell matching ####
umap_coproj = coembed@reductions$umap@cell.embeddings
ac_cells <- colnames(coembed)[coembed$tech == "ATAC"]
rna_cells <- colnames(coembed)[coembed$tech == "RNA"]
umap.rna = umap_coproj[rna_cells, ]
umap.atac = umap_coproj[ac_cells, ]
final_matching <- data.table(atac_cell = ac_cells)
final_matching$atac_cell <- as.character(final_matching$atac_cell)

dist0 <- pracma::distmat(umap.atac, umap.rna)
final_matching$rna_cell <- sapply(1:nrow(umap.atac), function(x) names(which.min(dist0[x, ])))


if(F){
  ## slower version -- not used
  final_matching$rna_cell<-sapply(final_matching$atac_cell, function(ac) {
    ac_umap <- umap_coproj[c(ac, rna_cells), c(1,2)]
    knn_k <- 1
    knn.res = FNN::get.knn(ac_umap, k = knn_k)
    knn.idx <- knn.res$nn.index[1,1]
    return(rownames(ac_umap)[knn.idx])
  })
  
}

saveRDS(final_matching, "Seurat_Objects/Integrated/atac_rna_cell_matching.rds")

## prepare metacell expression and peak accessiblity ####
## find nearest k = 10 cells for each cell in its orginal assay (RNA or ATAC)
K = 10
seurat.rna <- FindNeighbors(seurat.rna, k.param = K, reduction = 'pca')
knn.mat.rna = (seurat.rna@graphs$RNA_nn > 0)

seurat.atac <- FindNeighbors(seurat.atac, k.param = K, reduction = 'pca')
knn.mat.atac = (seurat.atac@graphs$ATAC_nn > 0)

all(rowSums(knn.mat.rna) == K)
all(rowSums(knn.mat.atac) == K)

smooth.rna = seurat.rna@assays$RNA@data %*% t(knn.mat.rna)
smooth.atac = seurat.atac@assays$ATAC@data %*% t(knn.mat.atac)


saveRDS(smooth.rna, file = "Seurat_Objects/Integrated/rna_metacell_expr.rds")
saveRDS(smooth.atac, file = "Seurat_Objects/Integrated/atac_metacell_access.rds")

## regression ####
coembed =readRDS(file = 'Seurat_Objects/Integrated/seurat_18MLLr_40Kcoembed.rds')
smooth.rna = readRDS("Seurat_Objects/Integrated/rna_metacell_expr.rds")
smooth.atac = readRDS("Seurat_Objects/Integrated/atac_metacell_access.rds")
smooth.rna = smooth.rna/10
smooth.atac = smooth.atac/10

final_matching = readRDS("Seurat_Objects/Integrated/atac_rna_cell_matching.rds")
smooth.rna = smooth.rna[, final_matching$rna_cell]
smooth.atac = smooth.atac[, final_matching$atac_cell]

## construct gene peak affinity binary matrix
gene_ann = fread('MetaData/gene_ann_hg38.txt')
gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]

final.peaks = rownames(smooth.atac)
final.peaks = sapply(final.peaks, function(x) unlist(strsplit(x, ','))[1])
names(final.peaks) = NULL
rownames(smooth.atac) = final.peaks
access.frac = rowMeans(smooth.atac > 0)
final.peaks = final.peaks[access.frac > 0.01]

## filter peaks that accessible in less than 10% of all cell type
seurat.atac = subset(coembed, tech == 'ATAC')

peaks.mean.ctype <- sapply(unique(seurat.atac$projCtype), function(x){
  rowMeans(seurat.atac@assays$ATAC@data[, seurat.atac$projCtype == x] > 0)
})

rmax = apply(peaks.mean.ctype, 1,  max)
summary(rmax)
final.peaks = names(which(rmax > 0.05))
final.peaks = lapply(final.peaks, function(x) unlist(strsplit(x, ','))[1])
final.peaks = do.call('c', final.peaks)

## focus on DEGs
degs1 = read.table("/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/PT_Ctype0/DEGs_betweenPTCtype0.txt")
degs1 = data.table(degs1)
degs1 = degs1[avg_logFC > 0.25 & p_val_adj < 0.05]
degs1 = unique(as.character(degs1$gene))

deg.dir = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/stagewise_DEG_5HDProjection/'
dfiles = dir(deg.dir)
dfiles = dfiles[grepl(dfiles, pattern = 'Alloutput')]
degs2 = lapply(dfiles, function(x){
  tmp = read.table(paste0(deg.dir, x))
  tmp = data.table(tmp, keep.rownames = T)
  tmp = tmp[abs(avg_logFC) > 0.25 & p_val_adj < 0.05]
})
degs2 = do.call('rbind', degs2)
degs2 = unique(degs2$rn)

degs3 = fread('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/DEG/Within_BTraj/DEGs_within_B_trajectory.txt')
degs3 = degs3[avg_logFC > 0.25 & p_val_adj < 0.05]
degs3 = unique(degs3$gene)

degs = unique(c(degs1, degs2, degs3))

## gene-peak-affinity map
gene2peak.map <- get_gene2peak_map(gene.list = data.table('gene_name' = degs),
                                   peak_names = final.peaks, 
                                   gene_ann = gene_ann,
                                   distal_dist = 5e5)

peaks.used = colnames(gene2peak.map)
degs = degs[degs %in% rownames(gene2peak.map)]

## do regression gene by gene
smooth.rna = smooth.rna[degs, ]
smooth.atac = smooth.atac[peaks.used, ]
smooth.rna = data.frame(as.matrix(smooth.rna))
#smooth.atac = data.frame(as.matrix(smooth.atac))


stime = Sys.time()
regr.list = list()
for(x in degs){
  exprs <- as.numeric(smooth.rna[x, ])
  names(exprs) <- NULL
  peaks = names(which(gene2peak.map[x, ] == 1))
  if(length(peaks) == 1) {
    covrs <- smooth.atac[peaks, ]
  }else{
    covrs <- t(smooth.atac[peaks, ])
  }
  rdata <- data.frame(cbind(exprs, covrs))
  res <- coef(summary(lm(exprs ~ ., data = rdata)))
  colnames(res)[4] <- 'P_value'
  regr.list[[x]] <- res
}
etime = Sys.time()
etime-stime
names(regr.list) = degs
saveRDS(regr.list, "EP_Prediction/regrRes4ep_prediction.rds")



## summarize/filter loops ####
regr.sum <- lapply(degs, function(t){
  x = regr.list[[t]]
  x = x[, c(1, 4)]
  x = data.frame(x)
  x = data.table(x, keep.rownames = T)
  x$gene_name = t
  return(x)
})

regr.sum = do.call('rbind', regr.sum)

regr.sum[, 'p_val_adj' := pmin(1, P_value*nrow(regr.sum))]
regr.sum$fdr = p.adjust(regr.sum$P_value, method = 'fdr')

regr.filtered = regr.sum[fdr < 0.05 & Estimate > 0.1 & grepl(rn, pattern = '^chr')]

regr.filtered$peak_name = sapply(regr.filtered$rn, function(x) gsub('.', '-', x, fixed = T)  )
regr.filtered$rn <- NULL

## to visualize on ucsc genome browser
## (promoter side: closest peak to TSS -- gene level)
gene_ann.deg = gene_ann[gene_name %in% regr.filtered$gene_name, ]

gene_ann.deg[, 'promoter_start' := Tss - 1000]
gene_ann.deg[, 'promoter_end' := Tss + 1000]

setkey(gene_ann.deg, gene_name)
regr.filtered[, 'start' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_start]
regr.filtered[, 'end' := gene_ann.deg[J(regr.filtered$gene_name)]$promoter_end]
regr.filtered[, 'chr' := gene_ann.deg[J(regr.filtered$gene_name)]$chr]
regr.filtered[, 'promoter_pos' := paste(chr, start, end, sep = '-')]

## filter otherend not overlapping with promoters
tss_ann = fread('MetaData/transcript_ann_hg38.txt')
tss_ann = tss_ann[gene_biotype %in% c('protein_coding', 'lincRNA', 'miRNA')]
tss_ann[, 'Tss' := ifelse(strand == '+', start, end)]

# any peak within promoter region
peak.ann = annPeak2Gene(peaks.used, tss_ann, 2000)
setkey(peak.ann, peak_name)
peaks.nprom = peak.ann[nchar(gene_name) == 0]$peak_name
peaks.prom = peak.ann[nchar(gene_name) > 0]$peak_name

regr.filtered.ep = regr.filtered[peak_name %in% peaks.nprom]

## assign nearest peak to promoter
gene_list <- subset(regr.filtered.ep, 
                   select = c(gene_name, chr, start, end)) %>%
           .[!duplicated(.)]

gene_list2peak = geneOverlapPeak(gene_list, peak_names = peaks.prom,
                                 mid_dist = 1000)
gene_list2peak = gene_list2peak[peak_name != 'Not_Found']
setkey(gene_list2peak, gene_name)

regr.filtered.ep = regr.filtered.ep[gene_name %in% gene_list2peak$gene_name]
regr.filtered.ep[, 'promoter_peak' := gene_list2peak[J(regr.filtered.ep$gene_name)]$peak_name]


regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, peak_name,
                                                       P_value, p_val_adj, fdr, Estimate))
names(regr.filtered.ep)[4] = 'enhancer_peak'

regr.filtered.ep[, 'chr1' := unlist(strsplit(promoter_peak, '-'))[1], by = promoter_peak]
regr.filtered.ep[, 'start1' := as.integer(unlist(strsplit(promoter_peak, '-'))[2]), by = promoter_peak]
regr.filtered.ep[, 'end1' := as.integer(unlist(strsplit(promoter_peak, '-'))[3]), by = promoter_peak]

regr.filtered.ep[, 'chr2' := unlist(strsplit(enhancer_peak, '-'))[1], by = enhancer_peak]
regr.filtered.ep[, 'start2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[2]), by = enhancer_peak]
regr.filtered.ep[, 'end2' := as.integer(unlist(strsplit(enhancer_peak, '-'))[3]), by = enhancer_peak]

regr.filtered.ep[, 'ep_dist' := abs(start1 + end1 - start2 - end2)/2]
regr.filtered.ep = subset(regr.filtered.ep, select = c(gene_name, promoter_pos, 
                                                       promoter_peak, enhancer_peak, ep_dist,
                                                       P_value, p_val_adj, fdr, Estimate))

fwrite(regr.filtered.ep, file = 'EP_Prediction/regrRes4_EP_overall.txt',
        sep = '\t')



