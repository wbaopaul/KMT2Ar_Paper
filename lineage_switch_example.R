source('scDataAnalysis_Utilities.R')

## 1. integration (and label transfer) -- atac ####
## intergrated seurat obtained by scATAC-pro (integrate module)
## then use script_transfer_label.R
seurat.obj <- readRDS('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/HTAN_1979_2741/integrated/seurat_obj_VFACS.rds')
DimPlot(seurat.obj, group.by = 'sample')
seurat.obj$sample[seurat.obj$sample == 'sample1'] = 'HTAN1979'
seurat.obj$sample[seurat.obj$sample == 'sample2'] = 'HTAN2524'
seurat.obj$sample[seurat.obj$sample == 'sample3'] = 'HTAN2578'
seurat.obj$sample[seurat.obj$sample == 'sample4'] = 'HTAN2741'
saveRDS(seurat.obj, file = 'Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_1979_All4TimePoints.rds')


seurat.obj2 <- readRDS('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/HTAN_2184_2263/integrated/seurat_obj_VFACS.rds')
DimPlot(seurat.obj2, group.by = 'sample')
seurat.obj2$cell_pre = ''
seurat.obj2$cell_pre[seurat.obj2$sample == 'sample1'] = 'CD33_Pos'
seurat.obj2$cell_pre[seurat.obj2$sample == 'sample2'] = 'CD33_Neg'

seurat.obj2$sample[seurat.obj2$sample == 'sample3'] = 'HTAN2184'
seurat.obj2$sample[seurat.obj2$sample %in% c('sample1', 'sample2')] = 'HTAN2263'

cnames = colnames(seurat.obj2)
cnames = paste0(seurat.obj2$cell_pre, cnames)
cnames = sapply(cnames, function(x) gsub('.1', '', x, fixed = T))
names(cnames) = NULL
seurat.obj2 = RenameCells(seurat.obj2, new.names = cnames)

saveRDS(seurat.obj2, file = 'Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_2184_2263.rds')


## 2. prepare for projection -- atac ####
## pool matrices for projection
## < 2184_2263 ####
## combine two sample
seurat1 <- readRDS('Seurat_Objects/scATAC/OtherLeuk/seurat_ALL2184_4projection.rds')
seurat2 <- readRDS('Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_2263_4projection.rds')
seurat1$sample = 'HTAN2184'
seurat2$sample = 'HTAN2263'

#rename cell
cnames1 = colnames(seurat1)
cnames2 = colnames(seurat2)
cnames2 = sapply(cnames2, function(x) gsub('Pos_', 'Pos', x))
cnames2 = sapply(cnames2, function(x) gsub('Neg_', 'Neg', x))
seurat2 = RenameCells(seurat2, new.names = cnames2)

seurat.comb = merge(seurat1, seurat2)

## added label 
seurat.obj <- readRDS('Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_2184_2263.rds')
mdata = subset(seurat.obj@meta.data, select = 'seurat_ctype')
seurat.comb = AddMetaData(seurat.comb, metadata = mdata)
saveRDS(seurat.comb, file = 'Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_2184_2263_4projection.rds')

## < 1979-2741 ####
mtx1 = read_mtx_scATACpro('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/OtherLeuk_Softlink/ALL1979/output/matrix4proj/matrix.mtx')
mtx2 = read_mtx_scATACpro('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/OtherLeuk_Softlink/HTAN2524/output/matrix4proj/matrix.mtx')
mtx3 = read_mtx_scATACpro('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/OtherLeuk_Softlink/HTAN2578/output/matrix4proj/matrix.mtx')
mtx4 = read_mtx_scATACpro('/mnt/isilon/tan_lab/yuw1/run_scATAC-pro/MLLr/OtherLeuk_Softlink/HTAN2741/output/matrix4proj/matrix.mtx')
shared.peaks = intersect(rownames(mtx1), rownames(mtx2))
shared.peaks = intersect(rownames(mtx3), shared.peaks)
shared.peaks = intersect(rownames(mtx4), shared.peaks)
mtx = cbind(mtx1[shared.peaks, ], mtx2[shared.peaks, ], 
            mtx3[shared.peaks, ], mtx4[shared.peaks, ])

samples = rep(c('HTAN1979', 'HTAN2524', 'HTAN2578', 'HTAN2741'),
              c(ncol(mtx1), ncol(mtx2), ncol(mtx3), ncol(mtx4)))

tss_ann = fread('/mnt/isilon/tan_lab/yuw1/scATAC-pro/annotation/hg38_tss.bed')
names(tss_ann)[c(1:4, 6, 7)] = c('chr', 'start', 'end', 'gene_name', 'strand',
                                 'gene_type')
mtx = assignGene2Peak(mtx, gene_ann = tss_ann)
seurat.comb = CreateSeuratObject(mtx, assay = 'ATAC')
seurat.comb$sample = samples

#add seurat ctype
seurat.obj <- readRDS('Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_1979_All4TimePoints.rds')
mdata = subset(seurat.obj@meta.data, select = 'seurat_ctype')
seurat.comb = AddMetaData(seurat.comb, metadata = mdata)


saveRDS(seurat.comb, file = 'Seurat_Objects/scATAC/OtherLeuk/seurat_HTAN_1979_All4TimePoints_4projection.rds')


## 3. projection -- atac ####
## < load hd data for projection ####
seuratAtacPath = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds'

seurat.hd <- readRDS(seuratAtacPath)

health_ann = seurat.hd@meta.data
health_ann = data.table(health_ann, keep.rownames = T)
setkey(health_ann, rn)
## revovery all the previous transformation done for thess healthy donor cells
vaps = VariableFeatures(seurat.hd)
feature.load = seurat.hd@reductions$pca@feature.loadings

## get the tf-idf 
mtx0 = 1*(seurat.hd@assays$ATAC@counts > 0)
mtx0 = mtx0[vaps, ]

npeaks <- colSums(x = mtx0)
tf <- t(t(mtx0) / npeaks)
idf <- ncol(mtx0) / rowSums(x = mtx0)
norm.mtx0 <- Diagonal(n = length(idf), x = idf) %*% tf
all(names(idf) == vaps)

load('Seurat_Objects/scATAC/coef_seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.RData')

npc = ncol(seurat.hd@reductions$pca@cell.embeddings)
npc = 20
set.seed(42) 
umap_model = uwot::umap(as.matrix(seurat.hd@reductions$pca@cell.embeddings[, 1:npc]),
                        n_neighbors = 30, min_dist = 0.3,
                        ret_model = T, metric = 'cosine')


## < proj ####
sampleNames = c('HTAN_1979_All4TimePoints', 'HTAN_2184_2263')

dir0 = 'Seurat_Objects/scATAC/OtherLeuk/'
cell.map = NULL
for(sampleName in sampleNames){
  seurat4ann = readRDS(paste0(dir0, '/seurat_', sampleName, '_4projection.rds'))
  mllr_ann1 = data.table(seurat4ann@meta.data, keep.rownames = T)
  
  mtx = seurat4ann@assays$ATAC@counts
  setkey(mllr_ann1, rn)
  
  ## scale the data using coefficent learned in healthy donor cells
  shared.peaks = intersect(rownames(mtx), vaps)
  coef.lm1 = coef.lm[shared.peaks, ]
  centers.vaps1 = centers.vaps[shared.peaks]
  sd.vaps1 = sd.vaps[shared.peaks]
  idf1 = idf[shared.peaks]
  
  mtx1 = mtx[shared.peaks, ]
  covMat = cbind(rep(1, ncol(mtx1)))
  mtx1 = t(t(mtx1)/colSums(mtx1))  # tf
  mtx1 = Diagonal(n = length(idf1), x = idf1) %*% mtx1
  mtx1.scale = (mtx1 - coef.lm1 %*% t(covMat) - centers.vaps1)/sd.vaps1
  
  
  U = feature.load[shared.peaks, 1:npc]
  
  pca1 = t(mtx1.scale) %*% U ##projected pcs
  
  ## to be consistent with seurat seed
  umap1 = uwot::umap_transform(as.matrix(pca1), umap_model)
  
  rownames(umap1) = rownames(pca1)
  umap1 = umap1[rownames(umap1) %in% mllr_ann1$rn, ]
  mllr_ann1$Ctype = mllr_ann1$seurat_ctype
  mllr_ann1[, 'Ctype' := ifelse(Ctype =='Blasts', paste0(sample, '_', 'Blasts'),
                                Ctype)]
  
  Ctype1 = mllr_ann1[rownames(umap1), ]$Ctype
  
  if(sampleName == 'HTAN_2184_2263'){
    myColors <- c("HTAN2184_Blasts" = "#fb6a4a",
                  "HTAN2263_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "T" = "#A45628")
  }else{
    myColors <- c("HTAN1979_Blasts" = "#fc9272",
                  "HTAN2524_Blasts" = "#fb6a4a",
                  "HTAN2578_Blasts" = "#ef3b2c",
                  "HTAN2741_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "Mature-B" = "#4EAE49",
                  "NK" = "#FF7E00",
                  "pDC" = "#989898",
                  "T" = "#A45628")
  }
  myColors['Healthy'] = '#cccccc'
  pdata = rbind(umap_model$embedding, umap1)
  pdata = data.table(pdata)
  pdata = cbind(pdata, c(rep('Healthy', nrow(umap_model$embedding)), 
                         Ctype1))
  names(pdata) = c('UMAP_1', 'UMAP_2', 'ctype')
  
  set.seed(2020)
  pdata$ctype = factor(pdata$ctype, levels = c('Healthy', sort(unique(Ctype1))) )
  pp1 <- ggplot(pdata[sample(1:nrow(pdata), 20000)], aes(x = UMAP_1, y = UMAP_2, col = ctype)) + 
    geom_point(size = 0.2) + 
    scale_colour_manual(name = "annotation", values = myColors) + theme_classic()
  pp1
  ggsave(pp1, filename = paste0('Figures/scATAC/proj2Healthy_final/proj_direct_', sampleName,
                                '.eps'), device = 'eps', width = 9, height = 7)
  
  
  ## assign the cell type from nearest healthy donor cells
  rownames(umap_model$embedding) = colnames(seurat.hd)
  ss.sub = pracma::distmat(umap1, umap_model$embedding)
  
  mllr2health <- sapply(1:nrow(umap1), function(x) names(which.min(ss.sub[x, ])))
  
  cell.map.tmp <- data.table('patient_bc' = rownames(umap1),
                             'healthy_bc' = mllr2health)
  cell.map.tmp[, 'ctype_healthy' := health_ann[J(cell.map.tmp$healthy_bc)]$seurat_ctype]
  cell.map.tmp[, 'sample_healthy' := health_ann[J(cell.map.tmp$healthy_bc)]$sample]
  cell.map.tmp[, 'ctype_patient' := mllr_ann1[J(cell.map.tmp$patient_bc)]$Ctype]
  cell.map.tmp[, 'sample_patient' := sampleName]
  cell.map = rbind(cell.map, cell.map.tmp)
  
}


map2fname = 'MetaData/scATAC/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
write.table(cell.map, file = map2fname, sep = '\t',
            row.names = F, quote = F)

## < summarize projection results -- atac ####
map2fname = 'MetaData/scATAC/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
cell.map = fread(map2fname)
cell.map = cell.map[grepl(ctype_patient, pattern = 'Blast')]
b.types = c('CLP', 'Pre-pro-B', 'Pro-B',
            'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B')
m.types = c('pDC', 'cDC', 'DC-Progenitor', 
            'GMP', 'Mono')
cell.map$lineage = 'Others'
cell.map[, 'lineage' := ifelse(ctype_healthy %in% m.types, 'M-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% b.types, 'B-lineage',
                               lineage)]

cell.map = cell.map[lineage != 'Others']
cell.map[, 'n' := .N, by = list(ctype_patient, lineage) ]
cell.map[, 'N' := .N, by = ctype_patient ]

cell.map.pdata = subset(cell.map, select = c('sample_patient', 'ctype_patient', 'lineage', 'n', 'N')) %>%
  .[!duplicated(.)]
cell.map.pdata[, 'frac' := round(n/N, 4)]
cell.map.pdata$perc = cell.map.pdata$frac * 100

myColors = c('M-lineage' = '#1F78B4',
             'B-lineage' = '#A6CEE3')
p1 <- ggplot(cell.map.pdata[sample_patient == "HTAN_2184_2263"],
       aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%'))) +
  theme_classic() + xlab('') + ylab('Percentage')
p2 <- ggplot(cell.map.pdata[sample_patient != "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%')),
            nudge_x = 0.1, nudge_y = -0.1) +
  theme_classic() + xlab('') + ylab('Percentage')

ggsave(p1, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN_2184_2263_atac.eps', 
       device = 'eps', width = 5, height = 5)
ggsave(p2, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN1979_All4TimePoints_atac.eps', 
       device = 'eps', width = 8, height = 5)



## < summarize projection results (updated with progenitor) -- atac ####
map2fname = 'MetaData/scATAC/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
cell.map = fread(map2fname)
cell.map = cell.map[grepl(ctype_patient, pattern = 'Blast')]
b.types = c('CLP', 'Pre-pro-B', 'Pro-B',
            'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B')
m.types = c('pDC', 'cDC', 'DC-Progenitor', 
            'GMP', 'Mono')
cell.map$lineage = 'Others'
cell.map[, 'lineage' := ifelse(ctype_healthy %in% m.types, 'M-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% b.types, 'B-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% c('LMPP'), 
                               'LMPP', lineage)]

cell.map = cell.map[lineage != 'Others']
cell.map[, 'n' := .N, by = list(ctype_patient, lineage) ]
cell.map[, 'N' := .N, by = ctype_patient ]

cell.map.pdata = subset(cell.map, select = c('sample_patient', 
                                             'ctype_patient', 'lineage', 'n', 'N')) %>%
  .[!duplicated(.)]
cell.map.pdata[, 'frac' := round(n/N, 4)]
cell.map.pdata$perc = cell.map.pdata$frac * 100

myColors = c('M-lineage' = '#1F78B4',
             'B-lineage' = '#A6CEE3',
             'LMPP' = "#33A02C")
cell.map.pdata$lineage = factor(cell.map.pdata$lineage, 
                                levels = c('B-lineage', 'M-lineage', 
                                            'LMPP'))

p1 <- ggplot(cell.map.pdata[sample_patient == "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%'))) +
  theme_classic() + xlab('') + ylab('Percentage')
p2 <- ggplot(cell.map.pdata[sample_patient != "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%')),
            nudge_x = 0.1, nudge_y = -0.1) +
  theme_classic() + xlab('') + ylab('Percentage')

ggsave(p1, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN_2184_2263_atac_updated.eps', 
       device = 'eps', width = 5, height = 5)
ggsave(p2, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN1979_All4TimePoints_atac_updated.eps', 
       device = 'eps', width = 8, height = 5)
saveRDS(cell.map.pdata, file = 'MetaData/M_lineage_LMPP_frac_lineageSwitch_atac.rds')

## 4. projection -- rna ####
seurat.rna <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')

## revovery all the previous transformation done for thess healthy donor cells
vegs = VariableFeatures(seurat.rna)
feature.load = seurat.rna@reductions$pca@feature.loadings
load('Seurat_Objects/scRNA/coef_hd.RData')
set.seed(42) 
npc = 20
umap_model = uwot::umap(as.matrix(seurat.rna@reductions$pca@cell.embeddings[, 1:npc]),
                        n_neighbors = 30, min_dist = 0.3,
                        ret_model = T, metric = 'cosine')
health_ann = seurat.rna@meta.data
health_ann = subset(health_ann, select = c('sample', 'nCount_RNA', 'nFeature_RNA', 'perc.mito',  'Ctype'))
health_ann = data.table(health_ann, keep.rownames = T)
setkey(health_ann, rn)

dir0 = 'Seurat_Objects/scRNA/OtherLeuk/'
cell.map = NULL
sampleNames = c('HTAN_1979_All4TimePoints', 'HTAN_2184_2263')

for(sampleName in sampleNames){
  seurat.mllr = readRDS(paste0(dir0, 'seurat_', sampleName, '_doubletRemoved.rds'))
  mllr_ann1 = seurat.mllr@meta.data
  mllr_ann1 = data.table(mllr_ann1, keep.rownames = T)
  setkey(mllr_ann1, rn)
  
  seurat.mllr = NormalizeData(seurat.mllr)
  
  ## scale the data using coefficent learned in healthy donor cells
  shared.genes = intersect(rownames(seurat.mllr), vegs)
  coef.lm1 = coef.lm[shared.genes, ]
  centers.vegs1 = centers.vegs[shared.genes]
  sd.vegs1 = sd.vegs[shared.genes]
  
  
  mtx1 = seurat.mllr@assays$RNA@data[shared.genes, ]
  mdata1 = seurat.mllr@meta.data
  perc.mito1 = mdata1[, 'perc.mito']
  nCount_RNA1 = mdata1[, 'nCount_RNA']
  covMat = cbind(rep(1, length(perc.mito1)), perc.mito1, nCount_RNA1)
  mtx1.scale = (mtx1 - coef.lm1 %*% t(covMat) - centers.vegs1)/sd.vegs1
  
  
  U = feature.load[shared.genes, 1:20]
  
  pca1 = t(mtx1.scale) %*% U ##projected pcs
  
  
  ## to be consistent with seurat seed
  umap1 = uwot::umap_transform(as.matrix(pca1), umap_model)
  
  rownames(umap1) = rownames(pca1)
  umap1 = umap1[rownames(umap1) %in% mllr_ann1$rn, ]
  
  mllr_ann1$Ctype = mllr_ann1$Ctype0
  mllr_ann1[, 'Ctype' := ifelse(Ctype =='Blasts', paste0(sample, '_', 'Blasts'),
                                Ctype)]
  mllr_ann1[, 'Ctype' := gsub('CD33_', '', Ctype)]
  mllr_ann1[, 'Ctype' := gsub('Neg_', '', Ctype)]
  
  Ctype1 = mllr_ann1[rownames(umap1), ]$Ctype
  
  if(sampleName == 'HTAN_2184_2263'){
    myColors <- c("HTAN2184_Blasts" = "#fb6a4a",
                  "HTAN2263_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "T" = "#A45628")
  }else{
    myColors <- c("HTAN1979_Blasts" = "#fc9272",
                  "HTAN2524_Blasts" = "#fb6a4a",
                  "HTAN2578_Blasts" = "#ef3b2c",
                  "HTAN2741_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "Mature-B" = "#4EAE49",
                  "NK" = "#FF7E00",
                  "pDC" = "#989898",
                  "T" = "#A45628")
  }
  myColors['Healthy'] = '#cccccc'
  
  
  
  pdata = rbind(umap_model$embedding, umap1)
  pdata = data.table(pdata)
  pdata = cbind(pdata, c(rep('Healthy', nrow(umap_model$embedding)), 
                         Ctype1))
  names(pdata) = c('UMAP_1', 'UMAP_2', 'ctype')
  
  
  set.seed(2020)
  
  pdata$ctype = factor(pdata$ctype, levels = c('Healthy', sort(unique(Ctype1))) )
  pp1 <- ggplot(pdata[sample(1:nrow(pdata), 20000)], aes(x = UMAP_1, y = UMAP_2, col = ctype)) + 
    geom_point(size = 0.2) + 
    scale_colour_manual(name = "annotation", values = myColors) +
    theme_classic()
  
  ggsave(pp1, filename = paste0('Figures/scRNA/proj2Healthy/final/proj_direct_', sampleName,
                                '.eps'), device = 'eps', width = 9, height = 7)
  
  
  ## assign the cell type from nearest healthy donor cells
  rownames(umap_model$embedding) = colnames(seurat.rna)
  ss.sub = pracma::distmat(umap1, umap_model$embedding)
  mllr2health <- sapply(1:nrow(umap1), function(x) names(which.min(ss.sub[x, ])))
  
  cell.map.tmp <- data.table('patient_bc' = rownames(umap1),
                             'healthy_bc' = mllr2health)
  
  cell.map.tmp[, 'ctype_healthy' := health_ann[J(cell.map.tmp$healthy_bc)]$Ctype]
  cell.map.tmp[, 'sample_healthy' := health_ann[J(cell.map.tmp$healthy_bc)]$sample]
  cell.map.tmp[, 'ctype_patient' := mllr_ann1[J(cell.map.tmp$patient_bc)]$Ctype]
  cell.map.tmp[, 'sample_patient' := sampleName]
  
  
  cell.map = rbind(cell.map, cell.map.tmp)
  
}
map2fname = 'MetaData/scRNA/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
write.table(cell.map, file = map2fname, sep = '\t',
            row.names = F, quote = F)

## < summarize projection results -- rna ####
map2fname = 'MetaData/scRNA/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
cell.map = fread(map2fname)
cell.map = cell.map[grepl(ctype_patient, pattern = 'Blast')]
b.types = c( 'CLP', 'Pre-pro-B', 'Pro-B',
            'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B')
m.types = c('pDC', 'cDC', 'DC-Progenitor', 
            'GMP', 'Mono')
cell.map$lineage = 'Others'
cell.map[, 'lineage' := ifelse(ctype_healthy %in% m.types, 'M-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% b.types, 'B-lineage',
                               lineage)]

cell.map = cell.map[lineage != 'Others']
cell.map[, 'n' := .N, by = list(ctype_patient, lineage) ]
cell.map[, 'N' := .N, by = ctype_patient ]

cell.map.pdata = subset(cell.map, select = c('sample_patient', 'ctype_patient', 'lineage', 'n', 'N')) %>%
  .[!duplicated(.)]
cell.map.pdata[, 'frac' := round(n/N, 4)]
cell.map.pdata$perc = cell.map.pdata$frac * 100

myColors = c('M-lineage' = '#1F78B4',
             'B-lineage' = '#A6CEE3')
p1 <- ggplot(cell.map.pdata[sample_patient == "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%'))) +
  theme_classic() + xlab('') + ylab('Percentage')
p2 <- ggplot(cell.map.pdata[sample_patient != "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%')),
            nudge_x = 0.1, nudge_y = -0.1) +
  theme_classic() + xlab('') + ylab('Percentage')

ggsave(p1, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN_2184_2263_rna.eps', 
       device = 'eps', width = 5, height = 5)
ggsave(p2, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN1979_All4TimePoints_rna.eps', 
       device = 'eps', width = 8, height = 5)



## < summarize projection results (add progenitor) -- rna ####
map2fname = 'MetaData/scRNA/TwoLineageSwitchExample_Cells_proj2healthyCells.csv'
cell.map = fread(map2fname)
cell.map = cell.map[grepl(ctype_patient, pattern = 'Blast')]
b.types = c( 'CLP', 'Pre-pro-B', 'Pro-B',
             'Pre-B', 'Immature-B', 'Mature-B', 'Plasma-B')
m.types = c('pDC', 'cDC', 'DC-Progenitor', 
            'GMP', 'Mono')
cell.map$lineage = 'Others'
cell.map[, 'lineage' := ifelse(ctype_healthy %in% m.types, 'M-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% b.types, 'B-lineage',
                               lineage)]
cell.map[, 'lineage' := ifelse(ctype_healthy %in% c('LMPP'), 
                               'LMPP', lineage)]

cell.map = cell.map[lineage != 'Others']
cell.map[, 'n' := .N, by = list(ctype_patient, lineage) ]
cell.map[, 'N' := .N, by = ctype_patient ]

cell.map.pdata = subset(cell.map, select = c('sample_patient', 
                                             'ctype_patient', 'lineage', 'n', 'N')) %>%
  .[!duplicated(.)]
cell.map.pdata[, 'frac' := round(n/N, 4)]
cell.map.pdata$perc = cell.map.pdata$frac * 100

myColors = c('M-lineage' = '#1F78B4',
             'B-lineage' = '#A6CEE3',
             'LMPP' = "#33A02C")
cell.map.pdata$lineage = factor(cell.map.pdata$lineage, 
                                levels = c('B-lineage', 'M-lineage', 
                                            'LMPP'))
p1 <- ggplot(cell.map.pdata[sample_patient == "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%'))) +
  theme_classic() + xlab('') + ylab('Percentage')
p2 <- ggplot(cell.map.pdata[sample_patient != "HTAN_2184_2263"],
             aes(x = ctype_patient, y = perc, fill = lineage)) + 
  geom_bar(stat = 'identity', position = 'dodge') + 
  scale_fill_manual(name = "annotation", values = myColors) +
  geom_text(aes(x = ctype_patient, y = perc, label = paste0(perc, '%')),
            nudge_x = 0.1, nudge_y = -0.1) +
  theme_classic() + xlab('') + ylab('Percentage')

ggsave(p1, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN_2184_2263_rna_updated.eps', 
       device = 'eps', width = 5, height = 5)
ggsave(p2, filename = 'Figures/summary/LineageSwitch/frac_lineage_HTAN1979_All4TimePoints_rna_updated.eps', 
       device = 'eps', width = 8, height = 5)

saveRDS(cell.map.pdata, file = 'MetaData/M_lineage_LMPP_frac_lineageSwitch_rna.rds')

## 5. plot signature genes ####
plotExprFrac <- function(seurat.obj, feature = 'MPO', assay = 'RNA',
                         group.name = 'sample',
                         plot.group.order){
  feature.exprs = seurat.obj[[assay]]@counts[feature, ]
  g.types = unique(seurat.obj@meta.data[, group.name])
  exprs.g <- lapply(g.types, function(x) feature.exprs[seurat.obj[[group.name]] == x] )
  names(exprs.g) = g.types
  
  dd.dt = data.table('group' = g.types,
                     'frac' = sapply(g.types, function(x) mean(exprs.g[[x]] > 0)))
  dd.dt$group = factor(dd.dt$group, levels = plot.group.order)
  dd.dt$frac = round(dd.dt$frac, 3)
  ymax = min(1, max(dd.dt$frac)+0.05)
  p0 <- ggplot(data = dd.dt, aes(x = group, y = frac, fill = group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab('Fraction of expression') + ylim(0, ymax) +
    geom_text(aes(label = frac, x = group, y = frac), 
              position = position_dodge(width = 0.8), vjust = -0.6) 
  return(p0)
  
}
sampleNames = c('HTAN_1979_All4TimePoints', 'HTAN_2184_2263')
sampleName = sampleNames[1]
seurat.mllr = readRDS(paste0('Seurat_Objects/scRNA/OtherLeuk/seurat_', sampleName, '_doubletRemoved.rds'))
  
seurat.mllr = subset(seurat.mllr, Ctype0 == 'Blasts')
seurat.mllr$Ctype = paste0(seurat.mllr$sample, '_Blasts')
seurat.mllr$Ctype = sapply(seurat.mllr$Ctype, function(x) gsub('CD33_|CD33_Neg_', '', x))
p0 <- StackedVlnPlot(seurat.mllr, features = c('CD19', 'HOXA9', 'MEIS1', "VPREB1", "IGLL1",
                                          "CD79A", 'DNTT', 'MPO', 'CEBPA', 'CD33'),
               group.by = 'Ctype', 
               myColors = brewer.pal(length(unique(seurat.mllr$Ctype)), 
                                     name = 'Paired')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.6))
ggsave(p0, filename = paste0('Figures/summary/LineageSwitch/vlnplot_sele_genes4Blasts_', sampleName, '.eps'), 
       device = 'eps', width = 5, height = 8)

plotExprFrac(seurat.mllr, feature = 'DNTT', assay = 'RNA',
                         group.name = 'Ctype',
                         plot.group.order = unique(seurat.mllr$Ctype))


## 6. update projection plot for atac, using uwot v0.1.5 for consistent umap plot ####
## < load hd data for projection ####
seuratAtacPath = 'Seurat_Objects/scATAC/seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.rds'

seurat.hd <- readRDS(seuratAtacPath)

health_ann = seurat.hd@meta.data
health_ann = data.table(health_ann, keep.rownames = T)
setkey(health_ann, rn)
## revovery all the previous transformation done for thess healthy donor cells
vaps = VariableFeatures(seurat.hd)
feature.load = seurat.hd@reductions$pca@feature.loadings

## get the tf-idf 
mtx0 = 1*(seurat.hd@assays$ATAC@counts > 0)
mtx0 = mtx0[vaps, ]

npeaks <- colSums(x = mtx0)
tf <- t(t(mtx0) / npeaks)
idf <- ncol(mtx0) / rowSums(x = mtx0)
norm.mtx0 <- Diagonal(n = length(idf), x = idf) %*% tf
all(names(idf) == vaps)

load('Seurat_Objects/scATAC/coef_seurat_9HDsamplesPlusThreeHD_TFIDF_usingDAPs1.RData')

npc = ncol(seurat.hd@reductions$pca@cell.embeddings)
npc = 20
set.seed(42) 
umap_model = uwot::umap(as.matrix(seurat.hd@reductions$pca@cell.embeddings[, 1:npc]),
                        n_neighbors = 30, min_dist = 0.3,
                        ret_model = T, metric = 'cosine')


## < proj ####
sampleNames = c('HTAN_1979_All4TimePoints', 'HTAN_2184_2263')

dir0 = 'Seurat_Objects/scATAC/OtherLeuk/'
cell.map = NULL
for(sampleName in sampleNames){
  seurat4ann = readRDS(paste0(dir0, '/seurat_', sampleName, '_4projection.rds'))
  mllr_ann1 = data.table(seurat4ann@meta.data, keep.rownames = T)
  
  mtx = seurat4ann@assays$ATAC@counts
  setkey(mllr_ann1, rn)
  
  ## scale the data using coefficent learned in healthy donor cells
  shared.peaks = intersect(rownames(mtx), vaps)
  coef.lm1 = coef.lm[shared.peaks, ]
  centers.vaps1 = centers.vaps[shared.peaks]
  sd.vaps1 = sd.vaps[shared.peaks]
  idf1 = idf[shared.peaks]
  
  mtx1 = mtx[shared.peaks, ]
  covMat = cbind(rep(1, ncol(mtx1)))
  mtx1 = t(t(mtx1)/colSums(mtx1))  # tf
  mtx1 = Diagonal(n = length(idf1), x = idf1) %*% mtx1
  mtx1.scale = (mtx1 - coef.lm1 %*% t(covMat) - centers.vaps1)/sd.vaps1
  
  
  U = feature.load[shared.peaks, 1:npc]
  
  pca1 = t(mtx1.scale) %*% U ##projected pcs
  
  ## to be consistent with seurat seed
  umap1 = uwot::umap_transform(as.matrix(pca1), umap_model)
  
  rownames(umap1) = rownames(pca1)
  umap1 = umap1[rownames(umap1) %in% mllr_ann1$rn, ]
  mllr_ann1$Ctype = mllr_ann1$seurat_ctype
  mllr_ann1[, 'Ctype' := ifelse(Ctype =='Blasts', paste0(sample, '_', 'Blasts'),
                                Ctype)]
  
  Ctype1 = mllr_ann1[rownames(umap1), ]$Ctype
  
  if(sampleName == 'HTAN_2184_2263'){
    myColors <- c("HTAN2184_Blasts" = "#fb6a4a",
                  "HTAN2263_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "T" = "#A45628")
  }else{
    myColors <- c("HTAN1979_Blasts" = "#fc9272",
                  "HTAN2524_Blasts" = "#fb6a4a",
                  "HTAN2578_Blasts" = "#ef3b2c",
                  "HTAN2741_Blasts" = "#cb181d",
                  "Monocyte" = "#367EB7",
                  "Mature-B" = "#4EAE49",
                  "NK" = "#FF7E00",
                  "pDC" = "#989898",
                  "T" = "#A45628")
  }
  myColors['Healthy'] = '#cccccc'
  pdata = rbind(umap_model$embedding, umap1)
  pdata = data.table(pdata)
  pdata = cbind(pdata, c(rep('Healthy', nrow(umap_model$embedding)), 
                         Ctype1))
  names(pdata) = c('UMAP_1', 'UMAP_2', 'ctype')
  
  set.seed(2020)
  pdata$ctype = factor(pdata$ctype, levels = c('Healthy', sort(unique(Ctype1))) )
  pp1 <- ggplot(pdata[sample(1:nrow(pdata), 20000)], aes(x = UMAP_1, y = UMAP_2, col = ctype)) + 
    geom_point(size = 0.2) + 
    scale_colour_manual(name = "annotation", values = myColors) + theme_classic()
  pp1
  ggsave(pp1, filename = paste0('Figures/scATAC/proj2Healthy_final_consistent/proj_direct_', sampleName,
                                '.eps'), device = 'eps', width = 9, height = 7)
  
  
  
}

