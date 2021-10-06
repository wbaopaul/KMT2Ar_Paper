source('scDataAnalysis_Utilities.R')

`%notin%` = Negate(`%in%`)

## load patient's normal data ####
seurat.atac = readRDS('Seurat_Objects/scATAC/seurat_pool_18MLLr_TFIDF_vap10000.rds')

seurat.atac = subset(seurat.atac, Ctype0 %notin% c('Blasts', 'Progenitors'))

seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_integrated_18MLLr_normal_final.rds')


## redo seurat of atac ####
seurat.atac <- doBasicSeurat_atac_updated(seurat.atac@assays$ATAC@counts,
                                          npc = 50, top.variable = 10000,
                                          norm_by = 'tf-idf', 
                                          vap.min.frac = 0.005,
                                          meta.data = seurat.atac@meta.data)

seurat.atac = RunUMAP(seurat.atac, dims = 1:50, verbose = F)
DimPlot(seurat.atac, group.by = 'Ctype0')


vaps0 = VariableFeatures(seurat.atac)
seurat.atac <- RunPCA(seurat.atac, npcs = 50, verbose = F, features = vaps0)
seurat.atac <- FindNeighbors(seurat.atac, dims = 1:50)
seurat.atac <- FindClusters(seurat.atac, resolusion = 0.2)

mtx = seurat.atac@assays$ATAC@counts
freq.pks = rowMeans(mtx > 0)
filtered.pks = names(which(freq.pks < 0.005))
## reselect
niter = 1
k = 0 
npc = 50
nvap = 5000

pp = list()
repeat{
  k = k + 1
  if(k > niter) break
  mtx_by_cls <- sapply(unique(seurat.atac$seurat_clusters), function(x) {
    cl_data <- 1* (mtx[, seurat.atac$seurat_clusters == x] > 0)
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
  seurat.atac[['ATAC']]@data = data0
  rm(data0, tmp, mtx.norm)
  
  VariableFeatures(seurat.atac) = sele.features
  
  seurat.atac <- ScaleData(seurat.atac, features = sele.features)
  seurat.atac <- RunPCA(seurat.atac, npc = npc, features = sele.features,
                        verbose = F)
  #seurat.atac <- regress_on_pca(seurat.atac)
  seurat.atac <- RunUMAP(seurat.atac, dims = 1:npc)
  
  pp[[k]] <- DimPlot(seurat.atac, group.by = 'sample', label = T)
  seurat.atac <- FindNeighbors(seurat.atac, reduction = 'pca', dims = 1:npc)
  seurat.atac <- FindClusters(seurat.atac, resolution = 0.2)
  
}


## transfer label ####
rn = rownames(mtx)
rownames(mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])

activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/yuw1/local_tools/annotation/GRCh38_genes.gtf',
                                            mtx, 
                                            include_body = T)
activity.matrix = activity.matrix[, colnames(activity.matrix) %in% colnames(seurat.atac)]
rm(mtx)


seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- FindVariableFeatures(seurat.atac)
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- ScaleData(seurat.atac)

DefaultAssay(seurat.atac) <- "ATAC"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'


## transfer label 

transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)


celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$Ctype_Final,
                                     weight.reduction = seurat.atac[["pca"]],
                                     dims = 1:ncol(seurat.atac[["pca"]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)
rm(transfer.anchors)
seurat.atac$seurat_ctype <- seurat.atac$predicted.id
seurat.atac$seurat_ctype_score_max <- seurat.atac$prediction.score.max

p1 <- DimPlot(seurat.atac, group.by = "predicted.id",
              label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells")
p2 <- DimPlot(seurat.rna, group.by = "Ctype_Final", label = TRUE,
              repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend() 
#p1
seurat.atac[["ACTIVITY"]] <- NULL ## don't save activity assay
saveRDS(seurat.atac, 'Seurat_Objects/scATAC/seurat_pool_18MLLr_normal.rds')

