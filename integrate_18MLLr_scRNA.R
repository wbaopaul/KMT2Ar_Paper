source('scDataAnalysis_Utilities.R')

## combine matrices ####



## load each individual scRNA seurat obj and 
## integrate by seurat3 ####
dir0 = 'Seurat_Objects/scRNA'

sampleNames = c('MLLr876533', 'MLLr882304', 'MLLr875706',
                'MLLr1154', 'MLLr870684', 'MLLr879583',
                'MLLr876545', 'MLLr874013', 'MLLr879440',
                'MLLr878289', 'MLLr871427', 'MLLr875703',
                'MLLr877476', 'MLLr881823', 'MLLr877780',
                'MLLr879339', 'MLLr878501', 'MLLr878516')
seurat.list = list()
for(id in sampleNamess){
  seurat.list[[sampleName]] <- readRDS(paste0(dir0, '/seurat_', sampleName.old,
                                              '_scRNA_doubletRemoved.rds'))
  seurat.list[[sampleName]]$sample = sampleName
}


seurat.obj <- integrateSeuratList(seurat.list, npc = 50, reg.var = NULL,
                                  reduction = 'pca', anchor.features = 3000)
rm(seurat.list)
seurat.obj <- RunUMAP(seurat.obj, dims = 1:50)
seurat.obj <- RunTSNE(seurat.obj, dims = 1:50)
seurat.obj <- FindNeighbors(seurat.obj, dims = 1:50, reduction = 'pca')
seurat.obj <- FindClusters(seurat.obj, res = 0.2)
DimPlot(seurat.obj, group.by = 'sample')

saveRDS(seurat.obj, 'Seurat_Objects/scRNA/seurat_integrated_18scRNA_VEG3000.rds')


## filtering variable genes #####
mtx = seurat.rna[['RNA']]@counts

freq_gene <- rowMeans(mtx > 0)
filter.genes = names(which(freq_gene < 0.002))
vegs = VariableFeatures(seurat.rna)

vegs.filtered = setdiff(vegs, filter.genes) 

seurat.rna <- RunPCA(seurat.rna, npcs = 50, verbose = F,
                      features = vegs.filtered)
seurat.rna <- RunUMAP(seurat.rna, dims = 1:50, verbose = F)

DimPlot(seurat.rna, group.by = 'sample', label = T)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:50, reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.2)


## reselect variable genes ####
niter = 2
k = 0
topn = 3000
npc = 50
repeat{
  k = k + 1
  if(k > niter) break
  clusters = as.character(seurat.rna$seurat_clusters)
  mtx_by_cls <- sapply(unique(clusters), function(x) {
    
    cl_data <- mtx[, clusters == x]
    
    Matrix::rowSums(cl_data)
    
  })
  mtx_by_cls.norm <- edgeR::cpm(mtx_by_cls, log = T, prior.count = 1)
  sds = sapply(1:nrow(mtx_by_cls.norm), function(x) sd(mtx_by_cls.norm[x, ]))
  names(sds) = rownames(mtx_by_cls.norm)
  sds = sort(sds, decreasing = T)
  sele.genes = names(sds[1:topn])
  sele.genes = setdiff(sele.genes, filter.genes)
  
  VariableFeatures(seurat.rna) <- sele.genes
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress =c('S.Score', 'G2M.Score', 'perc.mito',
                                             'HeatShock.Score1', 'nCount_RNA' ))
  seurat.rna <- RunPCA(seurat.rna, npcs = npc, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:npc)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:npc)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:npc, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.2)
  
}
DimPlot(seurat.rna, group.by = 'sample', label = T)

saveRDS(seurat.rna, file = 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds' )


## add manually annotated cell types information ####
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds' )
load('Input4VisCello/scRNA/relate2_seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated/State-2020-07-14_Ctype0.rda')
mdata = r_data$cmeta$df
mdata = subset(mdata, select = Ctype0)
seurat.rna <- AddMetaData(seurat.rna, metadata = mdata)
DimPlot(seurat.rna, group.by = 'Ctype0')
saveRDS(seurat.rna, file = 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds' )

