## transfer label from scRNA to scATAC using Seurat 

source('scDataAnalysis_Utilities.R')

## load scRNA data and cell type annotation ####
seurat.rna <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')
seuratAtacPath = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds'

seurat.atac = readRDS(seuratAtacPath)

getPalette = colorRampPalette(brewer.pal(9, "Paired"))
myColors = getPalette(17)
names(myColors) = c('cDC', 'CLP', 'DC_Progenitor', 'GMP',
                    'HSPC', 'Mature_B', 'LMPP', 'Immature_B',
                    'MEP', 'Mono', 'Pro-B', 'pDC', 'Plasma_B',
                    'Pre-B', 'T', 'NK', 'Pre-pro-B')

## use GAS = promote + gene body accessibility ####
atac.mtx = seurat.atac@assays$ATAC@counts
rn = rownames(atac.mtx)
rownames(atac.mtx) <- sapply(rn, function(x) unlist(strsplit(x, ','))[1])
activity.matrix = generate_gene_cisActivity('/mnt/isilon/tan_lab/yuw1/local_tools/annotation/GRCh38_genes.gtf',
                                            atac.mtx, 
                                            include_body = T)
activity.matrix = activity.matrix[, colnames(activity.matrix) %in% colnames(seurat.atac)]


#activity.matrix <- CreateGeneActivityMatrix(peak.matrix = atac.mtx, 
#                                            annotation.file = "/mnt/isilon/tan_lab/yuw1/local_tools/annotation/GRCh38_genes.gtf", 
#                                            seq.levels = c(1:22, "X", "Y"), upstream = 2000, verbose = TRUE)


seurat.atac[["ACTIVITY"]] <- CreateAssayObject(counts = activity.matrix)

DefaultAssay(seurat.atac) <- "ACTIVITY"
seurat.atac <- NormalizeData(seurat.atac)
seurat.atac <- ScaleData(seurat.atac, features = rownames(seurat.atac))
seurat.atac <- FindVariableFeatures(seurat.atac)

DefaultAssay(seurat.atac) <- "ATAC"
seurat.atac$tech = 'ATAC'
seurat.rna$tech = 'RNA'

## transfer label 
genes4anchors = VariableFeatures(object = seurat.rna)
#genes4anchors = NULL
transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.atac,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        query.assay = "ACTIVITY",
                                        reduction = "cca",
                                        k.anchor = 5)


celltype.predictions <- TransferData(anchorset = transfer.anchors,
                                     refdata = seurat.rna$Ctype,
                                     weight.reduction = seurat.atac[["pca"]],
                                     dims = 1:ncol(seurat.atac[["pca"]]),
                                     k.weight = 50)
celltype.predictions = subset(celltype.predictions, select = c('predicted.id', 'prediction.score.max'))
seurat.atac <- AddMetaData(seurat.atac, metadata = celltype.predictions)
#rm(transfer.anchors)

seurat.atac$seurat_ctype <- seurat.atac$predicted.id
seurat.atac$seurat_ctype_score_max <- seurat.atac$prediction.score.max

#seurat.atac.filtered <- subset(seurat.atac,
#                               subset = prediction.score.max > 0.3)


p1 <- DimPlot(seurat.atac, group.by = "predicted.id",
              label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
   scale_color_manual(values = myColors) 
p2 <- DimPlot(seurat.rna, group.by = "Ctype", label = TRUE,
              repel = TRUE) + ggtitle("scRNA-seq cells") + NoLegend() +
  scale_color_manual(values = myColors)
p1
#CombinePlots(plots = list(p1, p2))


seurat.atac[["ACTIVITY"]] <- NULL ## don't save activity assay
saveRDS(seurat.atac, seuratAtacPath)




