source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)


## using gene body + 2kb methylation matrix ####
## given a matrix list, cbind them using the union set of features, set NA for missing
cBind_union_features_NA <- function(mat_list){
  ff = rownames(mat_list[[1]])
  for(i in 2:length(mat_list)){
    ff = unique(union(ff, rownames(mat_list[[i]])))
  }
  ## make a mtx with full features
  mat_union = list()
  for(i in 1:length(mat_list)){
    mtx0 = mat_list[[i]]
    ff0 = setdiff(ff, rownames(mtx0))
    if(length(ff0) > 0 ) {
      tmp = as(matrix(NA, length(ff0), ncol(mtx0)), "sparseMatrix")
      rownames(tmp) = ff0
      tmp_mat = rbind(mtx0, tmp)
    }else{
      tmp_mat = mtx0
    }
    
    mat_union[[i]] = tmp_mat[order(rownames(tmp_mat)), ]
  }
  return(do.call('cbind', mat_union))
}
cBind_union_features_NA = cmpfun(cBind_union_features_NA)

impute_mtx <- function(mtx, impute.by = 'mean'){
  mtx1 = list()
  for(i in 1:ncol(mtx)){
    imtx = mtx[, i]
    if(impute.by == 'mean') imtx[is.na(imtx)] = mean(imtx[!is.na(imtx)])
    mtx1[[i]] = imtx
  }
  mtx1 = do.call('cbind', mtx1)
  colnames(mtx1) = colnames(mtx)
  return(mtx1)
}

dir0 = '/mnt/isilon/tan_lab/uzuny/projects/cptca/real_samples/data/snmc/LEUK/pilot//working_dir//Pilot//matrices/'
dir1 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_1/working_dir/'
dir2 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_2/working_dir/'
dir3 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/snmC/batch_3/working_dir/'

mtx.snmc = list()
mtx.snmc[['MLLr1154']] = as.matrix(readRDS(paste0(dir0, 'met_mat.genes_2kb.rds')))

for(idir in c(dir1, dir2, dir3)){
  sampleIDs = dir(idir)
  for(sampleID0 in sampleIDs){
    sampleName = paste0('MLLr', sampleID0)
    tmp = readRDS(paste0(idir, sampleID0, '/matrices/met_mat.genes_2kb.rds'))
    mtx.snmc[[sampleName]] <- tmp
  }
}
mtx.comb = cBind_union_features_NA(mtx.snmc)
samples = rep(names(mtx.snmc), sapply(mtx.snmc, ncol))


## filter by genes appeared in scRNA-seq data
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')

mtx.comb = mtx.comb[rownames(mtx.comb) %in% rownames(seurat.rna), ]

batches = rep('batch3', length(samples))
batches[samples == 'MLLr1154'] = 'Pilot'
batches[samples %in% c('MLLr877476', 'MLLr878516', 'MLLr882304')] = 'batch1'
batches[samples %in% c('MLLr874013', 'MLLr875706', 'MLLr879583')] = 'batch2'

mtx.imputed = impute_mtx(mtx.comb)
seurat.snmc = CreateSeuratObject(1 - log2(mtx.imputed+0.5), assay = 'snmc', 
                                 names.delim = '___')

## alternative formula
#seurat.snmc = CreateSeuratObject(log2(mtx.imputed+10^(-4))/log2(10^(-4)), assay = 'snmc', 
#                                 names.delim = '___')


seurat.snmc$batch = batches
seurat.snmc$sample = samples

seurat.snmc <- FindVariableFeatures(seurat.snmc)
seurat.snmc <- ScaleData(seurat.snmc, do.scale = T, do.center = T, 
                         vars.to.regress = 'nCount_snmc')
seurat.snmc <- RunPCA(seurat.snmc, npcs = 30, verbose = F)
seurat.snmc <- RunUMAP(seurat.snmc, dims = 1:30, verbose = F)
seurat.snmc <- FindNeighbors(seurat.snmc, dims = 1:30)
seurat.snmc <- FindClusters(seurat.snmc, resolution = 0.2)
p0 <- DimPlot(seurat.snmc, group.by = 'sample', label = T)
p1 <- DimPlot(seurat.snmc, group.by = 'batch', label = T)
FeaturePlot(seurat.snmc, features = 'nCount_snmc')

## < label transfer from scRNA ####
seurat.rna = subset(seurat.rna, Ctype0 == 'Blasts' & sample %in% unique(seurat.snmc$sample))

## transfer label 
transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.snmc,
                                        features = VariableFeatures(seurat.rna),
                                        reference.assay = "RNA",
                                        query.assay = "snmc",
                                        reduction = "cca",
                                        k.anchor = 5)
celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$Ctype_Stage,
                                     weight.reduction = seurat.snmc[["pca"]],
                                     dims = 1:ncol(seurat.snmc[["pca"]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, 
                              select = c('predicted.id', 'prediction.score.max'))
names(celltype.predictions) = c("seurat_ctype", "seurat_ctype_score_max")

seurat.snmc <- AddMetaData(seurat.snmc, metadata = celltype.predictions)

p1 <- DimPlot(seurat.snmc, group.by = "seurat_ctype",
              label = TRUE, repel = TRUE) + ggtitle("snmC-seq ") 
p2 <- DimPlot(seurat.rna, group.by = "Ctype_Stage", label = TRUE,
              repel = TRUE) + ggtitle("scRNA-seq ") + NoLegend() 
table(seurat.snmc$seurat_ctype)

saveRDS(seurat.snmc, file = 'Seurat_Objects/snmC/seurat_snmc_gene_2kb.rds')


## < label transfer from patient scATAC-seq ####
seurat.atac = readRDS('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap10000.rds')
seurat.atac = subset(seurat.atac, sample %in% unique(seurat.snmc$sample) & Ctype0 == 'Blasts')

## use GAS = promote + gene body accessibility 
atac.mtx = seurat.atac@assays$ATAC@counts
rn = rownames(atac.mtx)
rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/yuw1/local_tools/annotation/GRCh38_genes.gtf',
                                            atac.mtx, 
                                            include_body = T)
activity.matrix = activity.matrix[, colnames(activity.matrix) %in% colnames(seurat.atac)]
rm(atac.mtx)


seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)


DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- FindVariableFeatures(seurat.atac)
vegs = VariableFeatures(seurat.atac)
seurat.atac <- ScaleData(seurat.atac, features = vegs)

seurat.atac$tech = 'ATAC'
seurat.snmc$tech = 'snmC'

transfer.anchors <- FindTransferAnchors(reference = seurat.atac,
                                        query = seurat.snmc,
                                        features = vegs,
                                        reference.assay = "ACTIVITY",
                                        query.assay = "snmc",
                                        reduction = "cca",
                                        k.anchor = 5)
celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.atac$projCtype,
                                     weight.reduction = seurat.snmc[["pca"]],
                                     dims = 1:ncol(seurat.snmc[["pca"]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, 
                              select = c('predicted.id', 'prediction.score.max'))
names(celltype.predictions) = c("seurat_ctype_atac", "seurat_ctype_atac_score_max")

seurat.snmc <- AddMetaData(seurat.snmc, metadata = celltype.predictions)

p1 <- DimPlot(seurat.snmc, group.by = "seurat_ctype_atac",
              label = TRUE, repel = TRUE) + ggtitle("snmC-seq ") 
p2 <- DimPlot(seurat.atac, group.by = "projCtype", label = TRUE,
              repel = TRUE) + ggtitle("scATAC-seq ") + NoLegend() 
table(seurat.snmc$seurat_ctype_atac)

saveRDS(seurat.snmc, file = 'Seurat_Objects/snmC/seurat_snmc_gene_2kb.rds')



