source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)

## using gene body + 2kb as input ####

## < coembedding with scRAN -- not used ####
seurat.snmc <- readRDS('Seurat_Objects/snmC/seurat_snmc_gene_2kb.rds')
seurat.rna = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.rna = subset(seurat.rna, Ctype0 == 'Blasts' & sample %in% unique(seurat.snmc$sample))

seurat.snmc$tech = 'snmc'
seurat.rna$tech = 'RNA'

## transfer label 

genes4anchors = VariableFeatures(seurat.rna)
transfer.anchors <- FindTransferAnchors(reference = seurat.rna,
                                        query = seurat.snmc,
                                        features = genes4anchors,
                                        reference.assay = "RNA",
                                        query.assay = "snmc",
                                        reduction = "cca",
                                        k.anchor = 5)

#co-embedding
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.rna, assay = "RNA", slot = "data")[genes4anchors, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.snmc[["pca"]],
                           dims = 1:ncol(seurat.snmc[["pca"]]))

seurat.snmc[["RNA"]] <- imputation
coembed <- merge(x = seurat.rna, y = seurat.snmc)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
FeaturePlot(coembed, feature = 'pseudotime') 
saveRDS(coembed, file = 'Seurat_Objects/snmC/seurat_snmc_gene_body_2kb_coembedRNA.rds')

#1 to 1 matching
umap_coproj = coembed@reductions$umap@cell.embeddings
snmc_cells <- colnames(coembed)[coembed$tech == "snmc"]
rna_cells <- colnames(coembed)[coembed$tech == "RNA"]
umap.rna = umap_coproj[rna_cells, ]
umap.snmc = umap_coproj[snmc_cells, ]
final_matching <- data.table(snmc_cell = snmc_cells)
final_matching$snmc_cell <- as.character(final_matching$snmc_cell)

dist0 <- pracma::distmat(umap.snmc, umap.rna)
final_matching$rna_cell <- sapply(1:nrow(umap.snmc), 
                                  function(x) names(which.min(dist0[x, ])))

mdata.rna = data.table(seurat.rna@meta.data, keep.rownames = T)
mdata.snmc = data.table(seurat.snmc@meta.data, keep.rownames = T)
setkey(mdata.rna, rn)
setkey(mdata.snmc, rn)

final_matching[, 'projCtype' := mdata.rna[J(final_matching$rna_cell)]$projCtype]
final_matching[, 'pseudotime' := mdata.rna[J(final_matching$rna_cell)]$pseudotime]
final_matching[, 'sample' := mdata.snmc[J(final_matching$snmc_cell)]$sample]
final_matching[, 'batch' := mdata.snmc[J(final_matching$snmc_cell)]$batch]

saveRDS(final_matching, file = 'MetaData/snmC/cell_coembed_snmc_gene_body_2kb_referenceRNA.rds')


## < coembedding with scATAC ####
seurat.snmc <- readRDS('Seurat_Objects/snmC/seurat_snmc_gene_2kb.rds')
seurat.snmc <- readRDS('Seurat_Objects/snmC/seurat_snmc_gene_2kb_alterFormula.rds')
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
seurat.snmc$tech = 'snmc'

## transfer label 

genes4anchors = VariableFeatures(seurat.atac)
transfer.anchors <- FindTransferAnchors(reference = seurat.atac,
                                        query = seurat.snmc,
                                        features = genes4anchors,
                                        reference.assay = "ACTIVITY",
                                        query.assay = "snmc",
                                        reduction = "cca",
                                        k.anchor = 5)

#co-embedding
# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to

refdata <- GetAssayData(seurat.atac, assay = "ACTIVITY", slot = "data")[genes4anchors, ]

imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, 
                           weight.reduction = seurat.snmc[["pca"]],
                           dims = 1:ncol(seurat.snmc[["pca"]]))

seurat.snmc[["ACTIVITY"]] <- imputation
coembed <- merge(x = seurat.atac, y = seurat.snmc)

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- ScaleData(coembed, features = genes4anchors, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genes4anchors, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
FeaturePlot(coembed, feature = 'pseudotime') 
p1 <- DimPlot(coembed, group.by = 'tech')

ggsave(p1, filename = 'Figures/summary/umap_snmc_coembed_atac.eps', width = 7, height = 7,
       device = 'eps')

saveRDS(coembed, file = 'Seurat_Objects/snmC/seurat_snmc_gene_body_2kb_coembedATAC.rds')
#saveRDS(coembed, file = 'Seurat_Objects/snmC/seurat_snmc_gene_body_2kb_coembedATAC_alterFormula.rds')


#1 to 1 matching
umap_coproj = coembed@reductions$umap@cell.embeddings
snmc_cells <- colnames(coembed)[coembed$tech == "snmc"]
atac_cells <- colnames(coembed)[coembed$tech == "ATAC"]
umap.atac = umap_coproj[atac_cells, ]
umap.snmc = umap_coproj[snmc_cells, ]
final_matching <- data.table(snmc_cell = snmc_cells)
final_matching$snmc_cell <- as.character(final_matching$snmc_cell)

dist0 <- pracma::distmat(umap.snmc, umap.atac)
final_matching$atac_cell <- sapply(1:nrow(umap.snmc), 
                                   function(x) names(which.min(dist0[x, ])))

mdata.atac = data.table(seurat.atac@meta.data, keep.rownames = T)
mdata.snmc = data.table(seurat.snmc@meta.data, keep.rownames = T)
setkey(mdata.atac, rn)
setkey(mdata.snmc, rn)

final_matching[, 'projCtype' := mdata.atac[J(final_matching$atac_cell)]$projCtype]
final_matching[, 'pseudotime' := mdata.atac[J(final_matching$atac_cell)]$pseudotime]
final_matching[, 'sample' := mdata.snmc[J(final_matching$snmc_cell)]$sample]
final_matching[, 'batch' := mdata.snmc[J(final_matching$snmc_cell)]$batch]

saveRDS(final_matching, file = 'MetaData/snmC/cell_coembed_snmc_gene_body_2kb_referenceATAC.rds')
#saveRDS(final_matching, file = 'MetaData/snmC/cell_coembed_snmc_gene_body_2kb_referenceATAC_alterFormula.rds')


## add metadata infromation into seurat
final_matching = readRDS('MetaData/snmC/cell_coembed_snmc_gene_body_2kb_referenceATAC.rds')

getPalette = colorRampPalette(brewer.pal(9, "Paired"))
myColors = getPalette(17)
names(myColors) = paste0(c('cDC', 'CLP', 'DC-Progenitor', 'GMP',
                    'HSPC', 'Mature-B', 'LMPP', 'Immature-B',
                    'MEP', 'Mono', 'Pro-B', 'pDC', 'Plasma-B',
                    'Pre-B', 'T', 'NK', 'Pre-pro-B'), '-like')
b.types = paste0(c('CLP', 'Mature-B', 'Immature-B',
                   'Pro-B',  'Pre-pro-B', 'Pre-B'), '-like')

mdata = data.frame(final_matching)
rownames(mdata) = mdata$snmc_cell
seurat.snmc = AddMetaData(seurat.snmc, metadata = mdata)
p2 <- DimPlot(seurat.snmc, group.by = 'projCtype') + 
  scale_color_manual(values = myColors)

coembed$projCtype[final_matching$snmc_cell] = final_matching$projCtype

DimPlot(coembed, group.by = 'projCtype')
snmc.coembed <- subset(coembed, tech == 'snmc')

p3 <- DimPlot(snmc.coembed, group.by = 'projCtype') + 
  scale_color_manual(values = myColors)

ggsave(p2, filename = 'Figures/summary/umap_snmc_labelFromATAC.eps', width = 7, height = 7,
       device = 'eps')
ggsave(p3, filename = 'Figures/summary/umap_coembed_snmcOnly_labelFromATAC.eps', width = 7, height = 7,
       device = 'eps')

myColors['Others'] = '#cccccc'

seurat.snmc$projCtype_merge = 'Others'
seurat.snmc$projCtype_merge[seurat.snmc$projCtype %in% b.types] = seurat.snmc$projCtype[seurat.snmc$projCtype %in% b.types]
seurat.snmc$projCtype_merge = factor(seurat.snmc$projCtype_merge, 
                                     levels = c('CLP-like', 'Pre-pro-B-like', 
                                                'Pro-B-like', 'Pre-B-like',
                                                'Immature-B-like', 'Mature-B-like',
                                                'Others'))

snmc.coembed$projCtype_merge = 'Others'
snmc.coembed$projCtype_merge[snmc.coembed$projCtype %in% b.types] = snmc.coembed$projCtype[seurat.snmc$projCtype %in% b.types]
snmc.coembed$projCtype_merge = factor(snmc.coembed$projCtype_merge, 
                                     levels = c('CLP-like', 'Pre-pro-B-like', 
                                                'Pro-B-like', 'Pre-B-like',
                                                'Immature-B-like', 'Mature-B-like',
                                                'Others'))

p4 <- DimPlot(seurat.snmc, group.by = 'projCtype_merge') + 
  scale_color_manual(values = myColors)
p5 <- DimPlot(snmc.coembed, group.by = 'projCtype_merge') + 
  scale_color_manual(values = myColors)

ggsave(p4, filename = 'Figures/summary/umap_snmc_labelFromATAC.eps', 
       width = 7, height = 7, device = 'eps')
ggsave(p5, filename = 'Figures/summary/umap_coembed_snmcOnly_labelFromATAC.eps', 
       width = 7, height = 7, device = 'eps')



