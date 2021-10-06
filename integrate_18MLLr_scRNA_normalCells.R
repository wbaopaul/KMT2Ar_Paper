source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)

seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')

## for all normal cells ####
`%notin%` = Negate(`%in%`)
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna = subset(seurat.rna, Ctype0 != 'Blasts' & Ctype_updated %notin% c('HSPC1',
                                                                             'doublets'))
seurat.rna <- FindDoublets(seurat.rna, exp_rate = 0.02, PCs = 1:30)

seurat.rna = FindVariableFeatures(seurat.rna, nfeatures = 2000)

## < filter var genes ####
mtx = seurat.rna[['RNA']]@counts

freq_gene <- rowMeans(mtx > 0)
filter.genes = names(which(freq_gene < 0.005))
vegs = VariableFeatures(seurat.rna)

vegs.filtered = setdiff(vegs, filter.genes) 
VariableFeatures(seurat.rna) = vegs.filtered
seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                        vars.to.regress =c('perc.mito', 'nCount_RNA' ))
seurat.rna <- RunPCA(seurat.rna, npcs = 30, verbose = F,
                     features = vegs.filtered)
seurat.rna <- RunUMAP(seurat.rna, dims = 1:30, verbose = F)

DimPlot(seurat.rna, group.by = 'sample')

seurat.rna <- FindNeighbors(seurat.rna, dims = 1:30, reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.2)


## < reselect variable genes ####
niter = 2
k = 0
topn = 1000
npc = 30
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
                          vars.to.regress =c('perc.mito', 'nCount_RNA' ))
  seurat.rna <- RunPCA(seurat.rna, npcs = npc, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:npc)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:npc)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:npc, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.2)
  
}
DimPlot(seurat.rna, group.by = 'sample')

saveRDS(seurat.rna, file = 'Seurat_Objects/scRNA/seurat_pool_18MLLr_normal.rds')

## < using seurat integration --updated ####

seurat.list = SplitObject(seurat.rna, split.by = 'PresumedFusion')
norm.method = 'LogNormalize'

for (i in names(seurat.list)) {
  if(norm.method == 'SCT') seurat.list[[i]] <- SCTransform(seurat.list[[i]], verbose = FALSE)
  seurat.list[[i]] <- NormalizeData(seurat.list[[i]], verbose = FALSE)
}
mtx = seurat.rna@assays$RNA@counts
rs = rowSums(mtx > 0)

mllr.features <- SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 1500)
#seurat.list <- PrepSCTIntegration(object.list = seurat.list, anchor.features = mllr.features)

mllr.featuers = setdiff(mllr.features, names(which(rs < 100)))
# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(seurat.list) == "MLL-ENL")

mllr.anchors <- FindIntegrationAnchors(object.list = seurat.list, 
                                       normalization.method = norm.method, 
                                       anchor.features = mllr.features,
                                       reference = reference_dataset)
mllr.integrated <- IntegrateData(anchorset = mllr.anchors, 
                                 normalization.method = norm.method)
if(norm.method != 'SCT') mllr.integrated = ScaleData(mllr.integrated, 
                                                     vars.to.regress = c('perc.mito', 'nCount_RNA'))
mllr.integrated <- RunPCA(object = mllr.integrated, verbose = FALSE, npc = 30)
mllr.integrated <- RunUMAP(object = mllr.integrated, dims = 1:30)
DimPlot(mllr.integrated, group.by = 'sample')
DimPlot(mllr.integrated, group.by = 'Ctype_updated', label = T)
mllr.integrated <- FindNeighbors(mllr.integrated, reduction = 'pca', dims = 1:30)
mllr.integrated <- FindClusters(mllr.integrated, resolution = 0.1)

#DefaultAssay(mllr.integrated) <- "RNA"
FeaturePlot(mllr.integrated, features = c('CD14', 'FCGR3A','CD1C', 'CLEC4C'))
FeaturePlot(mllr.integrated, features = c('CD3G', 'CD4', 'CD8A', 'NCAM1'))

mllr.integrated <- FindClusters(mllr.integrated, resolution = 0.2)
mllr.integrated <- FindClusters(mllr.integrated, resolution = 0.4)
mllr.integrated <- FindClusters(mllr.integrated, resolution = 0.6)
mllr.integrated <- FindClusters(mllr.integrated, resolution = 0.8)
mllr.integrated <- FindClusters(mllr.integrated, resolution = 1)

saveRDS(mllr.integrated, 'Seurat_Objects/scRNA/seurat_integrated_18MLLr_normal_final.rds')


## update manually annotated cell type ####
mllr.integrated = readRDS('Seurat_Objects/scRNA/seurat_integrated_18MLLr_normal_final.rds')
load('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Input4VisCello/SaveStates/state_save_18MLLr_NormalCell_Final_Sept16')
ctype = subset(r_data$cmeta$df, select = 'Ctype_New')
mllr.integrated = AddMetaData(mllr.integrated, metadata = ctype, col.name = 'Ctype_Final')
mllr.integrated$Doublet_Singlet <- NULL

saveRDS(mllr.integrated, 'Seurat_Objects/scRNA/seurat_integrated_18MLLr_normal_final.rds')

## update manually annotated cell type for all 18 patinet MLLr ####
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
load('/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/Scripts/Input4VisCello/SaveStates/state_save_18MLLr_NormalCell_Final_Sept16')
ctype = subset(r_data$cmeta$df, select = 'Ctype_New')
seurat.rna$Ctype_Final = seurat.rna$Ctype0
seurat.rna$Ctype_Final[seurat.rna$Ctype_Final 
                       %notin% c('Blasts', 'Progenitors')] = 'doublets'
seurat.rna$Ctype_Final[rownames(ctype)] = ctype$Ctype_New
seurat.rna$Ctype_Final[seurat.rna$Ctype_Final == 'Progenitors'] = 'HSPC1'

seurat.rna$Ctype_Stage = seurat.rna$Ctype_Final
seurat.rna$Ctype_Stage[seurat.rna$Ctype_Stage == 'Blasts'] = seurat.rna$projCtype[seurat.rna$Ctype_Stage == 'Blasts']

saveRDS(seurat.rna, 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')

