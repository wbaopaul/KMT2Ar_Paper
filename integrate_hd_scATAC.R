source('scDataAnalysis_Utilities.R')

## global setting ####
getPalette = colorRampPalette(brewer.pal(9, "Paired"))
myColors = getPalette(17)
names(myColors) = c('cDC', 'CLP', 'DC_Progenitor', 'GMP',
                    'HSPC', 'Mature_B', 'LMPP', 'Immature_B',
                    'MEP', 'Mono', 'Pro-B', 'pDC', 'Plasma_B',
                    'Pre-B', 'T', 'NK', 'Pre-pro-B')
gene_ann = fread('GRCh38_gene_level.tsv')
gene_ann = gene_ann[!duplicated(gene_ann$gene_name)]
gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]

## select related feautures (from rn_names) to given genes
sele_peaks_gDEG <- function(degs_list, rn_names, 
                            gene_ann, distal_dist = 1e05 ){
  peak_names = sapply(rn_names, function(x) unlist(strsplit(x, ','))[1] )
  peaks <- tidyr::separate(data.table(peak_name = peak_names), 
                           col = peak_name, remove = F,
                           into = c('chr', 'start', 'end'))
  peaks$rn = rn_names
  class(peaks$start) = class(peaks$end) = 'integer'
  setkey(gene_ann, gene_name)
  setkey(peaks, rn)
  degs = data.table('rn' = degs_list)
  degs = degs[rn %in% gene_ann$gene_name]
  degs$chr = gene_ann[degs$rn]$chr
  degs$Tss = gene_ann[degs$rn]$Tss
  
  
  ## for each gene, get the corresponding peaks
  rname2gene <- lapply(1:nrow(degs), function(x) {
    chr0 = degs[x]$chr
    Tss0 = degs[x]$Tss
    
    peaks0 = peaks[chr == chr0]
    peaks0 = peaks0[abs(start/2 + end/2 - Tss0/2) <= distal_dist]
    return(peaks0$rn)
  } )
  
  return(rname2gene)
}


## 1. use TF-IDF -- add health donor samples from MPAL paper ####
## read seurat object outputted from scATAC-pro
ss = readRDS('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/output/HD_Uniq_Frag3000/integrated/seurat_obj_pool.rds') 

sampleNames = c("HD_7767_2111_CD34", "HD_7767_2111_LinNeg",
                "HD_7767_2111_live", "HD_7767_2153_CD34",
                "HD_7767_2153_lin_neg", "HD_7767_2153_live",
                "HD_7767_2689_CD34", "HD_7767_2689_lin_neg",
                "HD_7767_2689_live", 'Pub_HD1', 'Pub_HD2')
sample.match = data.table('sampleID' = paste0('sample', 1:11),
                          'sampleName' = sampleNames)
setkey(sample.match, sampleID)

ss$sample = sample.match[ss$sample]$sampleName


mtx = ss@assays$ATAC@counts
sampleIDs = ss$sample


rs = Matrix::rowSums(mtx > 0)
sele.pks1 = names(which(rs > 20))
mtx = mtx[sele.pks1, ]
batches = ifelse(grepl(sampleIDs, pattern = 'Pub_HD'), 'Pub_HD', 'HD')

freq_by_batch <- sapply(unique(batches), function(x){
  return(Matrix::rowMeans(mtx[, batches == x] > 0))
})

freq.max = rowMax(freq_by_batch)
freq.min = rowMin(freq_by_batch)
names(freq.max) = names(freq.min) =  rownames(freq_by_batch)
shared.peaks <- rownames(freq_by_batch)[freq.min > 0]
## save shared peaks
peaks2save = sapply(shared.peaks, function(x) unlist(strsplit(x, ','))[1])
peaks2save = tidyr::separate(data.table('V1' = peaks2save), col = V1, 
                             into = c('chr', 'start', 'end'))
write.table(peaks2save, file = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/output/HD_Uniq_Frag3000/peaks/shared_peaks_updated.bed',
            row.names = F, col.names = F, sep = '\t', quote = F)


# filter some features from vegs
rs = Matrix::rowSums(mtx > 0)
filter.pks1 = names(which(rs < (0.01 * ncol(mtx))))
filter.pks2 <- rownames(freq_by_batch)[freq.max > 4*freq.min]
filter.pks = unique(c(filter.pks1, filter.pks2))


mtx0 = mtx[shared.peaks, ]
seurat.atac0 = CreateSeuratObject(mtx0,
                                  assay = 'ATAC')
seurat.atac0$sample = sampleIDs
seurat.atac0$batch = batches

#seurat.atac0$nCount_ATAC = seurat.atac$nCount_ATAC
seurat.atac0@assays$ATAC@data = TF.IDF(mtx0)
seurat.atac0 <- FindVariableFeatures(seurat.atac0, 
                                    nfeatures = 5000)

#vegs = shared.peaks
vegs = VariableFeatures(seurat.atac0)
vegs = setdiff(vegs, filter.pks1)
vegs = setdiff(vegs, filter.pks2)

#manually add peaks associated with variable genes used in HD scRNA

VariableFeatures(seurat.atac0) <- vegs
seurat.atac0 <- ScaleData(seurat.atac0, 
                         features = VariableFeatures(seurat.atac0), 
                         vars.to.regress = c('nCount_ATAC'))
seurat.atac0 <- RunPCA(seurat.atac0, npc = 20, verbose = F, 
                      features = VariableFeatures(seurat.atac0))
seurat.atac0 <- RunUMAP(seurat.atac0, dims = 1:20, 
                       reduction = 'pca')
DimPlot(seurat.atac0, group.by = 'sample')

seurat.atac0 <- FindNeighbors(seurat.atac0, dims = 1:20, 
                             reduction = 'pca')
seurat.atac0 <- FindClusters(seurat.atac0, res = 0.4)
DimPlot(seurat.atac0, group.by = 'sample')

## reselect features 
## select variable features across clusters
library(edgeR)
niter = 1
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  clusters = as.character(seurat.atac0$seurat_clusters)
  mtx_by_cls <- sapply(unique(clusters), function(x) {
    
    cl_data <- mtx0[, clusters == x]
    
    Matrix::rowMeans(cl_data)
    
  })
  mtx_by_cls.norm <- cpm(mtx_by_cls, log = T, prior.count = 1)
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sele.features = names(which(sds >= sort(sds, decreasing = T)[2000]))
  sele.features <- setdiff(sele.features, filter.pks)
  mtx0.norm = TF.IDF(mtx0[sele.features, ])
  seurat.atac0@assays$ATAC@data[sele.features, ] <- mtx0.norm
  
  VariableFeatures(seurat.atac0) <- sele.features
  seurat.atac0 <- ScaleData(seurat.atac0, 
                           features = VariableFeatures(seurat.atac0), 
                           vars.to.regress = c('nCount_ATAC'))
  seurat.atac0 <- RunPCA(seurat.atac0, npc = 20, verbose = F, 
                        features = VariableFeatures(seurat.atac0))
  seurat.atac0 <- RunUMAP(seurat.atac0, dims = 1:20, 
                         reduction = 'pca')
  DimPlot(seurat.atac0, group.by = 'sample')
  seurat.atac0 <- FindNeighbors(seurat.atac0, dims = 1:20, 
                               reduction = 'pca')
  seurat.atac0 <- FindClusters(seurat.atac0, res = 0.4)
  
}
DimPlot(seurat.atac0, group.by = 'sample')


mtx <- seurat.atac@assays$ATAC@counts

freq_by_batch <- sapply(unique(seurat.atac$batch), function(x){
  return(Matrix::rowMeans(mtx[, seurat.atac$batch == x] > 0))
})

freq.max = rowMax(freq_by_batch)
freq.min = rowMin(freq_by_batch)
names(freq.max) = names(freq.min) =  rownames(freq_by_batch)
filter.pks2 <- rownames(freq_by_batch)[freq.max > 4*freq.min]
rs = Matrix::rowSums(mtx > 0)
filter.pks1 <- names(which(rs <= 0.01 * ncol(mtx)))

filter.pks = unique(union(filter.pks1, filter.pks2))

seurat.atac@assays$ATAC@data = TF.IDF(mtx)
seurat.atac <- FindVariableFeatures(seurat.atac, 
                                     nfeatures = 5000)


vegs = VariableFeatures(seurat.atac)
vegs = setdiff(vegs, filter.pks)

VariableFeatures(seurat.atac) <- vegs
seurat.atac <- ScaleData(seurat.atac, 
                          features = VariableFeatures(seurat.atac), 
                          vars.to.regress = c('nCount_ATAC'))
seurat.atac <- RunPCA(seurat.atac, npc = 20, verbose = F, 
                       features = VariableFeatures(seurat.atac))
seurat.atac <- RunUMAP(seurat.atac, dims = 1:20, 
                        reduction = 'pca')
DimPlot(seurat.atac, group.by = 'sample')

seurat.atac <- FindNeighbors(seurat.atac, dims = 1:20, 
                              reduction = 'pca')
seurat.atac <- FindClusters(seurat.atac, res = 0.4)


## reselect features 
## select variable features across clusters
library(edgeR)
niter = 2
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  clusters = as.character(seurat.atac$seurat_clusters)
  mtx_by_cls <- sapply(unique(clusters), function(x) {
    
    cl_data <- mtx[, clusters == x]
    
    Matrix::rowMeans(cl_data)
    
  })
  mtx_by_cls.norm <- cpm(mtx_by_cls, log = T, prior.count = 1)
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sele.features = names(which(sds >= sort(sds, decreasing = T)[2000]))
  sele.features <- setdiff(sele.features, filter.pks)
   
  mtx.norm = TF.IDF(mtx[sele.features, ])
  seurat.atac@assays$ATAC@data[sele.features, ] <- mtx.norm
  
  VariableFeatures(seurat.atac) <- sele.features
  seurat.atac <- ScaleData(seurat.atac, 
                            features = VariableFeatures(seurat.atac), 
                            vars.to.regress = c('nCount_ATAC'))
  seurat.atac <- RunPCA(seurat.atac, npc = 20, verbose = F, 
                         features = VariableFeatures(seurat.atac))
 
  seurat.atac <- RunUMAP(seurat.atac, dims = 1:20, 
                          reduction = 'pca')
  DimPlot(seurat.atac, group.by = 'sample')
  seurat.atac <- FindNeighbors(seurat.atac, dims = 1:20, 
                                reduction = 'pca')
  seurat.atac <- FindClusters(seurat.atac, res = 0.4)
  
}
DimPlot(seurat.atac, group.by = 'sample')


seuratFile = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusMPALHD_TFIDF_uniq_frag3000_RemovedDoublets_updated.rds'
saveRDS(seurat.atac, file = seuratFile )









## 4. add a new pub data (Pub_HD3) ####
pub3 = read_mtx_scATACpro('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/Pub_HD3/output/filtered_matrix/COMBINED/FILTER/matrix.mtx')
tss_ann = fread('/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_tss.bed')
names(tss_ann)[c(1:4, 6)] <- c('chr', 'start', 'end', 'gene_name',
                               'strand')
pub3 = assignGene2Peak(pub3, gene_ann = tss_ann)

mtx0 = seurat.atac@assays$ATAC@counts
shared.peaks = intersect(rownames(mtx0), rownames(pub3))

colnames(pub3) = paste0('Pub3_', colnames(pub3))
mtx = cbind(mtx0[shared.peaks, ], pub3[shared.peaks, ])
mdata = data.frame('sample' = c(seurat.atac$sample, 
                                rep('Pub_HD3', ncol(pub3))),
                   'batch' = c(seurat.atac$batch, 
                               rep('Pub_HDnew', ncol(pub3))),
                   stringsAsFactors = F)
rownames(mdata) = colnames(mtx)
rm(seurat.atac, pub3, mtx0)
seurat.atac.new <- CreateSeuratObject(mtx, assay = 'ATAC',
                                      meta.data = mdata)

seurat.atac.new$batch[seurat.atac.new$batch == 'Pub_HD'] = 'HD_MPAL'

freq_by_batch <- sapply(unique(seurat.atac.new$batch), function(x){
  return(Matrix::rowMeans(mtx[, seurat.atac.new$batch == x] > 0))
})

freq.max = apply(freq_by_batch, 1, max)
freq.min = apply(freq_by_batch, 1, min)
names(freq.max) = names(freq.min) =  rownames(freq_by_batch)
rs = rowMeans(mtx > 0)
filter.pks1 = names(which(rs < 0.01))
#filter.pks2 <- rownames(freq_by_batch)[freq.max > 4*freq.min]
filter.pks2 <- names(which(freq.max < 0.05))
filter.pks = unique(c(filter.pks1, filter.pks2))

seurat.atac.new[['ATAC']]@data = TF.IDF(seurat.atac.new[['ATAC']]@counts)

seurat.atac.new <- FindVariableFeatures(seurat.atac.new, 
                                    nfeatures = 8000)

vegs = VariableFeatures(seurat.atac.new)
#vegs = union(vegs, vfeatures)
vegs = setdiff(vegs, filter.pks)
VariableFeatures(seurat.atac.new) <- vegs

seurat.atac.new <- ScaleData(seurat.atac.new, 
                         features = VariableFeatures(seurat.atac.new), 
                         vars.to.regress = NULL)
seurat.atac.new <- RunPCA(seurat.atac.new, npc = 30, verbose = F, 
                      features = VariableFeatures(seurat.atac.new))
seurat.atac.new <- RunUMAP(seurat.atac.new, dims = 1:30, 
                       reduction = 'pca')
DimPlot(seurat.atac.new, group.by = 'sample')

seurat.atac.new <- FindNeighbors(seurat.atac.new, dims = 1:30, 
                             reduction = 'pca')
seurat.atac.new <- FindClusters(seurat.atac.new, res = 0.6)

library(edgeR)
niter = 1
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  clusters = as.character(seurat.atac.new$seurat_clusters)
  mtx_by_cls <- sapply(unique(clusters), function(x) {
    
    cl_data <- mtx[, clusters == x]
    
    Matrix::rowMeans(cl_data)
    
  })
  mtx_by_cls.norm <- cpm(mtx_by_cls, log = T, prior.count = 1)
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sele.features = names(which(sds >= sort(sds, decreasing = T)[2500]))
  sele.features <- setdiff(sele.features, filter.pks)
  mtx.norm = TF.IDF(mtx[sele.features, ])
  seurat.atac.new@assays$ATAC@data[sele.features, ] <- mtx.norm
  
  VariableFeatures(seurat.atac.new) <- sele.features
  seurat.atac.new <- ScaleData(seurat.atac.new, 
                           features = VariableFeatures(seurat.atac.new), 
                           vars.to.regress = NULL)
  seurat.atac.new <- RunPCA(seurat.atac.new, npc = 30, verbose = F, 
                        features = VariableFeatures(seurat.atac.new))
  
  seurat.atac.new <- RunUMAP(seurat.atac.new, dims = 1:30, 
                         reduction = 'pca')
  DimPlot(seurat.atac.new, group.by = 'sample')
  seurat.atac.new <- FindNeighbors(seurat.atac.new, dims = 1:30, 
                               reduction = 'pca')
  seurat.atac.new <- FindClusters(seurat.atac.new, res = 0.4)
  
}
DimPlot(seurat.atac.new, group.by = 'sample')

seurat.atac.new = readRDS(seuratFile)
mtx = seurat.atac.new@assays$ATAC@counts

freq_by_batch <- sapply(unique(seurat.atac.new$batch), function(x){
  return(Matrix::rowMeans(mtx[, seurat.atac.new$batch == x] > 0))
})

freq.max = apply(freq_by_batch, 1, max)
freq.min = apply(freq_by_batch, 1, min)
names(freq.max) = names(freq.min) =  rownames(freq_by_batch)
rs = rowMeans(mtx > 0)
filter.pks1 <- rownames(freq_by_batch)[freq.max > 5*freq.min]
filter.pks2 <- names(which(freq.max < 0.05))
filter.pks = unique(c(filter.pks1, filter.pks2))

seurat.atac.new[['ATAC']]@data = TF.IDF(seurat.atac.new[['ATAC']]@counts)
seurat.atac.new <- ScaleData(seurat.atac.new, 
                             features = rownames(seurat.atac.new), 
                             vars.to.regress = NULL)

seurat.atac.new <- FindVariableFeatures(seurat.atac.new, 
                                        nfeatures = 10000)

vegs = VariableFeatures(seurat.atac.new)
vegs = setdiff(vegs, filter.pks)
VariableFeatures(seurat.atac.new) <- vegs

seurat.atac.new <- RunPCA(seurat.atac.new, npc = 50, verbose = F, 
                          features = VariableFeatures(seurat.atac.new))
seurat.atac.new <- RunUMAP(seurat.atac.new, dims = 1:50, 
                           reduction = 'pca')
DimPlot(seurat.atac.new, group.by = 'sample')

seurat.atac.new <- FindNeighbors(seurat.atac.new, dims = 1:50, 
                                 reduction = 'pca')
seurat.atac.new <- FindClusters(seurat.atac.new, res = 0.6)

library(edgeR)
niter = 1

for(nvap in c(2000)){
  k = 0
  repeat{
    k = k + 1
    if(k > niter) break
    clusters = as.character(seurat.atac.new$seurat_clusters)
    mtx_by_cls <- sapply(unique(clusters), function(x) {
      
      cl_data <- mtx[, clusters == x]
      
      Matrix::rowMeans(cl_data)
      
    })
    mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
    sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
    names(sds) = rownames(mtx_by_cls.norm)
    sele.features = names(which(sds >= sort(sds, decreasing = T)[nvap]))
    sele.features <- setdiff(sele.features, filter.pks)
    mtx.norm = TF.IDF(mtx[sele.features, ])
    seurat.atac.new@assays$ATAC@data[sele.features, ] <- mtx.norm
    
    VariableFeatures(seurat.atac.new) <- sele.features
    seurat.atac.new <- RunPCA(seurat.atac.new, npc = 50, verbose = F, 
                              features = VariableFeatures(seurat.atac.new))
    
    seurat.atac.new <- RunUMAP(seurat.atac.new, dims = 1:50, 
                               reduction = 'pca')
    DimPlot(seurat.atac.new, group.by = 'sample')
    #seurat.atac.new <- FindNeighbors(seurat.atac.new, dims = 1:30, 
    #                                 reduction = 'pca')
    #seurat.atac.new <- FindClusters(seurat.atac.new, res = 0.4)
    
  }
  p1 <- DimPlot(seurat.atac.new, group.by = 'sample')
  
  tmpFile = paste0('Figures/test_umap_hd_atac_nvap', nvap, '.png')
  ggsave(p1, filename = tmpFile, device = 'png', width = 8, height = 8)
  
}

#nvap=2000 best; choose nvap=2000

seuratFile = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_vap2000.rds'
saveRDS(seurat.atac.new, file = seuratFile )


##  using daps as variable features ####
seuratFile = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_vap2000.rds'
seurat.atac.new = readRDS(seuratFile)
sele.features0 = VariableFeatures(seurat.atac.new)
all.daps <- FindAllMarkers(seurat.atac.new, test.use = 'wilcox',
                           logfc.threshold = 0, slot = 'data',
                           max.cells.per.ident = 200,
                           min.pct = 0.1, only.pos = T,
                           min.diff.pct = 0.05)

all.daps = data.table(all.daps, keep.rownames = T)
saveRDS(all.daps, file = 'Seurat_Objects/scATAC/DAPs_9HDsamplesPlusThreeHD_TFIDF_vap2000.rds')

all.daps = readRDS('Seurat_Objects/scATAC/DAPs_9HDsamplesPlusThreeHD_TFIDF_vap2000.rds')


mtx = seurat.atac.new@assays$ATAC@counts

freq_by_batch <- sapply(unique(seurat.atac.new$batch), function(x){
  return(Matrix::rowMeans(mtx[, seurat.atac.new$batch == x] > 0))
})

freq.max = apply(freq_by_batch, 1, max)
freq.min = apply(freq_by_batch, 1, min)
names(freq.max) = names(freq.min) =  rownames(freq_by_batch)
rs = rowMeans(mtx > 0)
filter.pks1 <- rownames(freq_by_batch)[freq.max > 5*freq.min]
filter.pks2 <- names(which(freq.max < 0.02))
filter.pks = unique(c(filter.pks1, filter.pks2))

sele.daps = all.daps[p_val_adj < 0.05 & avg_logFC > 0.1]
sele.features <- setdiff(sele.daps$gene, filter.pks)


mtx.norm0 <- TF.IDF(mtx[sele.features, ])
tmp <- mtx[setdiff(rownames(mtx), sele.features), ]
data0 <- rbind(mtx.norm0, tmp)
seurat.atac[['ATAC']]@data = data0[rownames(mtx), ]
rm(data0, tmp, mtx.norm0)

VariableFeatures(seurat.atac.new) <- sele.features
seurat.atac.new <- ScaleData(seurat.atac.new, features = sele.features)

seurat.atac.new <- RunPCA(seurat.atac.new, npc = 30, verbose = F, 
                          features = VariableFeatures(seurat.atac.new))

seurat.atac.new <- RunUMAP(seurat.atac.new, dims = 1:20, 
                           reduction = 'pca')
DimPlot(seurat.atac.new, group.by = 'sample')


seuratFile = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds'
saveRDS(seurat.atac.new, seuratFile)


