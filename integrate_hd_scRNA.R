## process data from healty donors 

source('scDataAnalysis_Utilities.R')


## use lognorm + gini for all sample ####
files = dir('Seurat_Objects/scRNA/HealthyDonor/')
files = files[grepl(files, pattern = 'seurat_Normal|seurat_HD')]
tmp <- readRDS(paste0('Seurat_Objects/scRNA/HealthyDonor/', files[1]))
tmp <- subset(tmp, nCount_RNA < 40000)
mtx <- tmp@assays$RNA@counts
samples = tmp$sample
perc.mito = tmp$perc.mito

for(file0 in files[-1]){
  tmp <- readRDS(paste0('Seurat_Objects/scRNA/HealthyDonor/', file0))
  tmp <- subset(tmp, nCount_RNA < 40000 )
  mtx0 <- tmp@assays$RNA@counts
  samples0 = tmp$sample
  perc.mito0 = tmp$perc.mito
  
  samples = c(samples, samples0)
  perc.mito = c(perc.mito, perc.mito0)
  
  mtx = t(Matrix.utils::rBind.fill(t(mtx), t(mtx0), fill = 0))
}

rm(tmp, mtx0)

seurat.rna <- CreateSeuratObject(mtx)
seurat.rna <- NormalizeData(seurat.rna)
seurat.rna$sample = samples
seurat.rna$perc.mito = perc.mito
seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 1000)

seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                        vars.to.regress = c('perc.mito', 'nCount_RNA'))
seurat.rna <- RunPCA(seurat.rna, npc = 20, verbose = F, 
                     features = VariableFeatures(seurat.rna))

seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
#seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, 
                            reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.2)


niter = 2
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  gini_genes <- ifg_select(seurat.rna@assays$RNA@counts, 
                           seurat.rna$seurat_clusters, 
                           gini_cut_qt = 0.9)$include_g
  
  VariableFeatures(seurat.rna) <- gini_genes
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress = c('perc.mito', 'nCount_RNA'))
  seurat.rna <- RunPCA(seurat.rna, npc = 20, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.1)
  
}
seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)


## try to remove potential doublets
seurat.rna <- FindDoublets(seurat.rna, PCs = 1:20, exp_rate = 0.01)

## Downsample progenitors ####
seurat.rna = subset(seurat.rna, Doublet_Singlet == 'Singlet')
seurat.rna$childID = sapply(seurat.rna$sample, 
                            function(x) unlist(strsplit(x, '_'))[2])


set.seed(2020)
ndown = 10000
id_live = grep(seurat.rna$sample, pattern = 'Live')
id_pg = setdiff(1:length(seurat.rna$sample), id_live)
id_sub = sample(id_pg, ndown)

##rm cd3 & CD19 coexpress cells
mtx = seurat.rna@assays$RNA@counts

sele_id = sort(c(id_live, id_sub))

mtx = seurat.rna@assays$RNA@counts[, sele_id]
samples = seurat.rna$sample[sele_id]
perc.mito = seurat.rna$perc.mito[sele_id]
sorts = seurat.rna$sort[sele_id]

seurat.rna <- CreateSeuratObject(mtx)
seurat.rna <- NormalizeData(seurat.rna)
seurat.rna$sample = samples
seurat.rna$perc.mito = perc.mito
seurat.rna$sort = sorts


seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 1000)

seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                        vars.to.regress = c('perc.mito', 'nCount_RNA'))
seurat.rna <- RunPCA(seurat.rna, npc = 20, verbose = F, 
                     features = VariableFeatures(seurat.rna))

seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
#seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, 
                            reduction = 'pca')
seurat.rna <- FindClusters(seurat.rna, res = 0.2)


niter = 2
k = 0
repeat{
  k = k + 1
  if(k > niter) break
  gini_genes <- ifg_select(seurat.rna@assays$RNA@counts, 
                           seurat.rna$seurat_clusters, 
                           gini_cut_qt = 0.9)$include_g
  
  VariableFeatures(seurat.rna) <- gini_genes
  
  seurat.rna <- ScaleData(seurat.rna, features = VariableFeatures(seurat.rna), 
                          vars.to.regress = c('perc.mito', 'nCount_RNA'))
  seurat.rna <- RunPCA(seurat.rna, npc = 20, verbose = F, 
                       features = VariableFeatures(seurat.rna))
  
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:20)
  #seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)
  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:20, reduction = 'pca')
  seurat.rna <- FindClusters(seurat.rna, res = 0.1)
  
}
seurat.rna <- RunTSNE(seurat.rna, dims = 1:20)

old.seurat <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')
mdata = subset(old.seurat@meta.data, select = Ctype)
seurat.rna <- AddMetaData(seurat.rna, metadata = mdata, 
                          col.name = 'Old_Ctype' )
seurat.rna$Old_Ctype = as.character(seurat.rna$Old_Ctype)
seurat.rna$Old_Ctype[is.na(seurat.rna$Old_Ctype)] = 'NA'

saveRDS(seurat.rna, 'Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')



