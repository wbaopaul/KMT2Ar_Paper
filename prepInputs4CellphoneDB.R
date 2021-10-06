library(data.table)
library(Seurat)
library(magrittr)
`%notin%` = Negate(`%in%`)

## < load seurat ####
seurat.mllr = readRDS('Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds')
seurat.mllr$sample = as.character(seurat.mllr$sample)
mllr_ann = data.table(seurat.mllr@meta.data, keep.rownames = T, stringsAsFactors = F)

sampleNames = unique(mllr_ann$sample)

table(mllr_ann[Ctype0 == 'Progenitors']$sample)
progenitor.g <- c('MLLr879440', 'MLLr875703', 'MLLr871427',
                   'MLLr1154', 'MLLr876545', 'MLLr870684', 
                  'MLLr875706', 'MLLr874013', 'MLLr882304')
younger.g <-  c('MLLr876533', 'MLLr882304', 'MLLr875706',
                'MLLr1154', 'MLLr870684', 'MLLr879583',
                'MLLr876545', 'MLLr874013', 'MLLr879440',
                'MLLr878289', 'MLLr871427')

## < pool input data into 4 groups ####
mllr_ann = subset(mllr_ann, Ctype_Stage %notin% c('T-like', 'NK-like', 'MEP-like', 'doublets',
                                           'Unknown'))
mllr_ann$Ctype_Stage[mllr_ann$Ctype_Stage %in% c( 'pDC-like', 'cDC-like', 'DC-Progenitor-like')] = 'DC-like'
mllr_ann.sele = mllr_ann


mllr_ann.prog = mllr_ann.sele[sample %in% progenitor.g]
mllr_ann.nprog = mllr_ann.sele[sample %notin% progenitor.g]
mllr_ann.younger = mllr_ann.sele[sample %in% younger.g]
mllr_ann.older = mllr_ann.sele[sample %notin% younger.g]


mllr_ann.prog = subset(mllr_ann.prog, select = c('rn', 'Ctype_Stage'))
mllr_ann.nprog = subset(mllr_ann.nprog, select = c('rn', 'Ctype_Stage'))
mllr_ann.younger = subset(mllr_ann.younger, select = c('rn', 'Ctype_Stage'))
mllr_ann.older = subset(mllr_ann.older, select = c('rn', 'Ctype_Stage'))
mllr_ann.sele = subset(mllr_ann.sele, select = c('rn', 'Ctype_Stage'))



names(mllr_ann.younger) = names(mllr_ann.prog) = names(mllr_ann.older) = c('Cell', 'cell_type')

names(mllr_ann.nprog) = names(mllr_ann.sele) = c('Cell', 'cell_type')

counts = seurat.mllr[['RNA']]@data[, mllr_ann.sele$Cell]
counts.prog = seurat.mllr[['RNA']]@data[, mllr_ann.prog$Cell]
counts.nprog = seurat.mllr[['RNA']]@data[, mllr_ann.nprog$Cell]
counts.younger = seurat.mllr[['RNA']]@data[, mllr_ann.younger$Cell]
counts.older = seurat.mllr[['RNA']]@data[, mllr_ann.older$Cell]


ntotal = 30000

## << With HSPC1 ####
set.seed(2019)
ids = sample((1:nrow(mllr_ann.prog)), ntotal)
mllr_ann.prog = mllr_ann.prog[ids, ]
ctype.rm.prog = names(which(table(mllr_ann.prog$cell_type) < 20))
mllr_ann.prog = mllr_ann.prog[cell_type %notin% ctype.rm.prog]
counts.prog = counts[, mllr_ann.prog$Cell]

dir.create('CellPhoneDB_18MLLr/PoolPatientsWHSPC1/', 
           showWarnings = F, recursive = T)
mfile.prog = paste0('CellPhoneDB_18MLLr/PoolPatientsWHSPC1/mdata.txt')
write.table(mllr_ann.prog, file = mfile.prog,
            sep = '\t', row.names = F, quote = F)

## combine blast into one group
mllr_ann.prog <- fread(mfile.prog)
mllr_ann.prog[grep(cell_type, pattern = 'like')]$cell_type = 'Blasts'
mfile.prog.blast = paste0('CellPhoneDB_18MLLr/PoolPatientsWHSPC1/mdata_blasts.txt')
write.table(mllr_ann.prog, file = mfile.prog.blast,
            sep = '\t', row.names = F, quote = F)


counts.prog = as.matrix(counts.prog)
counts.prog = data.table(counts.prog, keep.rownames = T)
names(counts.prog)[1] = 'Gene'
cfile.prog = paste0('CellPhoneDB_18MLLr/PoolPatientsWHSPC1/data.txt')
fwrite(counts.prog, file = cfile.prog,
       sep = '\t',  quote = F)
