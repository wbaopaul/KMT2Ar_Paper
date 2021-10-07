source('scDataAnalysis_Utilities.R')

hg38.ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/scMethod/GRCh38_gene_level.tsv')
cBind_union_features <- function(mtx_list){
  ff0 = ff = rownames(mtx_list[[1]])
  for(i in 2:length(mtx_list)){
    ff = unique(union(ff, rownames(mtx_list[[i]])))
  }
  if(all(ff0 == ff)) return(do.call('cbind', mtx_list))
  ## make a mtx with full features
  mat_union = list()
  for(i in 1:length(mtx_list)){
    mtx0 = mtx_list[[i]]
    ff0 = setdiff(ff, rownames(mtx0))
    if(length(ff0) > 0 ) {
      tmp = as(matrix(0, length(ff0), ncol(mtx0)), "sparseMatrix")
      rownames(tmp) = ff0
      mtx0 = rbind(mtx0, tmp)
    }
    mat_union[[i]] = mtx0[order(rownames(mtx0)), ]
  }
  return(do.call('cbind', mat_union))
}

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
myColors = getPalette(18)
names(myColors) =  c("PAYWJZ", "PAZGKI", "PAYUZM", "1154",  
                     "PAYKGI", "PAZBSZ", "PAYWKL", "PAYSBA", 
                      "PAZBLA", "PAYZLC", "PAYLNH", "PAYUZJ", 
                      "PAYYBG", "PAZFPH", "PAYYNY", "PAZBGV",
                      "PAYZVY", "PAYZWN ")


## analyzed pooled 18 MLLr data ####
## read and pool all mtx
dir0 = '/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/'
samples = list.dirs(dir0, recursive = F)
samples <- sapply(samples, basename ) %>% 
  sapply(., function(x) ifelse(grepl(x, pattern = 'MLL'), x, ''))
samples = samples[nchar(samples) > 0]
names(samples) = NULL
samples = samples[samples != 'MLLr_old_cells']

mtx_list = list()
sampleIDs = NULL
for(sample0 in samples){
  mtx_dir0 = paste0(dir0, sample0, 
                    '/output/filtered_matrix/COMBINED/FILTER/reConstruct_matrix/matrix.mtx')
  tmp = read_mtx_scATACpro(mtx_dir0)
  colnames(tmp) = paste0( sample0, colnames(tmp))
  mtx_list[[sample0]] = tmp
  
}

mtx = cBind_union_features(mtx_list)
rm(tmp, mtx_list)
rs = Matrix::rowMeans(mtx > 0)
mtx = mtx[rs > 0.001, ]

sampleIDs = sapply(colnames(mtx), function(x) 
  substr(x, 1, nchar(x) - 16))

sampleIDs = sapply(sampleIDs, function(x) gsub('MLL_', 'MLLr', x))
sampleIDs = sapply(sampleIDs, function(x) gsub('_scATAC', '', x))


tss_ann = fread('/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_tss.bed')
names(tss_ann)[c(1:4, 6, 7)] = c('chr', 'start', 'end', 'gene_name', 'strand',
                                 'gene_type')
tss_ann = tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding')] 
mtx = assignGene2Peak(mtx, gene_ann = tss_ann)

freq_by_sample <- sapply(unique(sampleIDs), function(x){
  return(Matrix::rowMeans(mtx[, sampleIDs == x] > 0))
})
freq.max = apply(freq_by_sample, 1, max)
freq.min = apply(freq_by_sample, 1, min)
names(freq.max) = names(freq.min) =  rownames(freq_by_sample)
filter.pks1 = names(which(freq.max < 0.01))

mtx = mtx[!rownames(mtx) %in% filter.pks1, ]  ## roughly filter
npc = 50
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.01))
nvap = 10000
seurat.mllr = doBasicSeurat_atac_updated(mtx, npc = npc, 
                                         norm_by = 'tf-idf',
                                         top.variable = nvap,
                                         regressOnPca = FALSE,
                                         reg.var = NULL, 
                                         excludePks.fromVAP = filtered.pks)
seurat.mllr$sample = sampleIDs

vaps1 = VariableFeatures(seurat.mllr)

seurat.mllr <- RunUMAP(seurat.mllr, dims = 1:npc)


DimPlot(seurat.mllr, group.by = 'sample', label = T) + 
  scale_color_manual(values = myColors)

## reselect variable features across clusters
seurat.mllr <- FindNeighbors(seurat.mllr, reduction = 'pca', 
                             dims = 1:npc)
seurat.mllr <- FindClusters(seurat.mllr, resolution = 0.4)
saveRDS(seurat.mllr, 
        file = 'Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF.rds')


niter = 2
k = 0 
repeat{
  k = k + 1
  if(k > niter) break
  mtx_by_cls <- sapply(unique(seurat.mllr$seurat_clusters), function(x) {
                         cl_data <- 1* (mtx[, seurat.mllr$seurat_clusters == x] > 0)
                         return(Matrix::rowSums(cl_data))
                 })
  
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.features = names(sds[1:nvap])
  sele.features = setdiff(sele.features, filtered.pks)
  mtx.norm = TF.IDF(mtx[sele.features, ])
  seurat.mllr[['ATAC']]@data[sele.features, ] <- mtx.norm
  VariableFeatures(seurat.mllr) = sele.features
  seurat.mllr <- ScaleData(seurat.mllr, features = sele.features)
  seurat.mllr <- RunPCA(seurat.mllr, npc = npc, features = sele.features)
  #seurat.mllr <- regress_on_pca(seurat.mllr)
  seurat.mllr <- RunUMAP(seurat.mllr, dims = 1:npc)
 
  
  seurat.mllr <- FindNeighbors(seurat.mllr, reduction = 'pca', dims = 1:npc)
  seurat.mllr <- FindClusters(seurat.mllr, resolution = 0.4)
  
}
DimPlot(seurat.mllr, group.by = 'sample', label = T) + 
  scale_color_manual(values = myColors)


saveRDS(seurat.mllr, 
        file = 'Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_old.rds')



## filtered variable peaks more stringent & using diff vaps ####
tmp = readRDS(file = 'Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_old.rds')
mtx = tmp[['ATAC']]@counts
rs = Matrix::rowMeans(mtx > 0)
filtered.pks = names(which(rs < 0.02))
npc = 50
nvap = 20000
seurat.mllr = doBasicSeurat_atac_updated(mtx,
                                         npc = npc, 
                                         top.variable = nvap,
                                         regressOnPca = FALSE,
                                         reg.var = NULL,
                                         meta.data = tmp@meta.data,
                                         excludePks.fromVAP = filtered.pks)
rm(tmp)
seurat.mllr <- RunUMAP(seurat.mllr, dims = 1:npc)
DimPlot(seurat.mllr, group.by = 'sample', label = T) +
  scale_color_manual(values = myColors)

seurat.mllr <- FindNeighbors(seurat.mllr, reduction = 'pca', dims = 1:npc)
seurat.mllr <- FindClusters(seurat.mllr, resolution = 0.2)
vaps0 = VariableFeatures(seurat.mllr)

## reselect
niter = 2
k = 0 

pp = list()
repeat{
  k = k + 1
  if(k > niter) break
  mtx_by_cls <- sapply(unique(seurat.mllr$seurat_clusters), function(x) {
    cl_data <- 1* (mtx[, seurat.mllr$seurat_clusters == x] > 0)
    return(Matrix::rowSums(cl_data))
  })
  
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.features = names(sds[1:nvap])
  sele.features = setdiff(sele.features, filtered.pks)
  
 
  mtx.norm = TF.IDF(mtx[sele.features, ])
  tmp <- mtx[setdiff(rownames(mtx), sele.features), ]
  data0 <- rbind(mtx.norm, tmp)
  seurat.mllr[['ATAC']]@data = data0
  rm(data0, tmp, mtx.norm)
  
  VariableFeatures(seurat.mllr) = sele.features
  
  seurat.mllr <- ScaleData(seurat.mllr, features = sele.features)
  seurat.mllr <- RunPCA(seurat.mllr, npc = npc, features = sele.features,
                        verbose = F)
  #seurat.mllr <- regress_on_pca(seurat.mllr)
  seurat.mllr <- RunUMAP(seurat.mllr, dims = 1:npc)
  
  pp[[k]] <- DimPlot(seurat.mllr, group.by = 'sample', label = T) + 
    scale_color_manual(values = myColors)
  seurat.mllr <- FindNeighbors(seurat.mllr, reduction = 'pca', dims = 1:npc)
  seurat.mllr <- FindClusters(seurat.mllr, resolution = 0.2)
  
}


fileName = paste0('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap', nvap, '.rds')
saveRDS(seurat.mllr, file = fileName)

## add ctype annotation as meta data
nvap = 10000
fileName = paste0('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap', nvap, '.rds')
seurat.atac = readRDS(fileName)
load('Input4VisCello/scATAC/relate2_seurat_pool_18MLLr_TFIDF_vap10000/State-2020-07-16_Ctype0.rda')
ctype0 = r_data$cmeta$df
ctype0 = subset(ctype0, select = 'Ctype0')
seurat.atac <- AddMetaData(seurat.atac, metadata = ctype0)
saveRDS(seurat.atac, file = fileName)
