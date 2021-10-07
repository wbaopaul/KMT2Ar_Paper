## ---------------------------------------------------------------##
## *** summarize longread seq final result **** ##
## ---------------------------------------------------------------##

source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`)

Add_fusionInf2Seurat <- function(sampleName, seurat.rna, nmismatch = 0,
                                 dtype = 'both'){
  sampleName1 = sampleName
  if(sampleName %in% c('PAZBLA', 'PAYLNH', 'PAYUZJ')) 
    sampleName1 = paste0(sampleName, '_Redo')
  
  if(sampleName %in% c('HD2111_CD34', 'HD2689_Live'))
    sampleName1 = unlist(strsplit(sampleName, '_'))[1]
  
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/LongRead/Final_Results/'
  read_file = paste0(dir0, sampleName1, '/', sampleName1, 
                     '.fastq.gz.demultiplexing.PASS.reads.txt')
  
  message(paste0(sampleName, ':'))
  
  ## read all long reads for used bc
  long.res = fread(read_file)
  names(long.res)[c(1:3, 7)] = c('readname', 'bc', 'mismatch', 'sec_mismatch')
  long.res = long.res[mismatch <= nmismatch]
  long.res = subset(long.res, select = c('readname', 'bc', 'mismatch', 'sec_mismatch'))
  #long.res = long.res[sec_mismatch > 4]
  n0 = nrow(long.res)
  
  ## filtered out reads that are not map2 KMT2A PCR region
  if(sampleName == 'PAYKGI') return(seurat.rna)
  reads.pcr = fread(paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_reads2PCR.txt'))
  long.res = long.res[readname %in% reads.pcr$readname, ]
  message(paste0(round(nrow(long.res)/n0 * 100, 1), '% reads in PCR!'))
  
 
  long.res$bc = paste0(sampleName, '_', long.res$bc)
  
  
  long.res[, 'nLongread' := .N, by = bc]
  
  
  ## note the bc with longread data
  tmp_long = subset(long.res, select = c(nLongread, bc)) %>% .[!duplicated(.)]
  setkey(tmp_long, bc)
  tmp_mdata = data.table('bc' = rownames(seurat.rna@meta.data[seurat.rna$sample == sampleName, ]))
  tmp_mdata[, 'nLongread' := tmp_long[J(tmp_mdata$bc)]$nLongread]
  tmp_mdata$nLongread[is.na(tmp_mdata$nLongread)] = 0
  
  seurat.rna$nLongread[tmp_mdata$bc] = tmp_mdata$nLongread
  
  long.res = long.res[bc %in% colnames(seurat.rna)]
  
  
  ## read fusion reads from alignqc
  if(dtype != 'longGF'){
    fusion.file.alignqc = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
    if(!file.exists(fusion.file.alignqc)) return(seurat.rna)
    fusion_reads.alignqc =  unique(fread(fusion.file.alignqc)$readname)
    long.fusion.alignqc = long.res[readname %in% fusion_reads.alignqc & bc %in% colnames(seurat.rna)]
    
    message(paste('AlignQC detected', nrow(long.fusion.alignqc), 'fusion reads in',
                  length(unique(long.fusion.alignqc$bc)), 'cells!'))
    
  }
  
 
  ## read fusion reads from longGF
  if(dtype != 'alignqc'){
    fusion.files = dir(paste0(dir0, sampleName1))
    fusion.file.sele = fusion.files[grepl(fusion.files, pattern = 'KMT2A|AP001267')]
    fusion.file.sele = fusion.file.sele[grepl(fusion.file.sele, pattern = 'list$')]
    fusion.file.sele = fusion.file.sele[grepl(fusion.file.sele, pattern = 'filtered.readname', fixed = T)]
    fusion_reads.longgf = lapply(fusion.file.sele, function(x) fread(paste0(dir0, sampleName1, '/', 
                                                                            x), header = F)$V1)
    fusion_reads.longgf = unique(do.call('c', fusion_reads.longgf))
    long.fusion.longgf = long.res[readname %in% fusion_reads.longgf & bc %in% colnames(seurat.rna)]
    
    message(paste('longGF detected', nrow(long.fusion.longgf), 'fusion reads in',
                  length(unique(long.fusion.longgf$bc)), 'cells!'))
  
  }
  
  if(dtype == 'alignqc') long.fusion = long.fusion.alignqc
  if(dtype == 'longGF')  long.fusion = long.fusion.longgf
  
  if(dtype == 'both') {
    long.fusion = long.fusion.alignqc[readname %in% long.fusion.longgf$readname]
    message(paste('Detected', nrow(long.fusion), 'fusion reads in',
                  length(unique(long.fusion$bc)), 'cells by both methods!'))
    
  }
  
  message('-----------------------------------')
  if(nrow(long.fusion) == 0) return(seurat.rna)
  
  ## write final read inf
  mdata0 = data.table(seurat.rna@meta.data, keep.rownames = T)
  mdata0$Ctype0[mdata0$Ctype0 %in% c('T', 'NKT')] = 'T/NK'
  setkey(mdata0, rn)
  long.fusion[, 'Ctype0' := mdata0[J(long.fusion$bc)]$Ctype0]
  write.table(subset(long.fusion[nLongread >= 1], select = c('readname', 'bc', 'Ctype0')),
              file = paste0('MetaData/scRNA/LongRead/Final_Fusions/', sampleName, '.txt'),
              row.names = F, quote = F, sep = '\t')
  
  long.res.sum = subset(long.fusion, select = c(bc))
  long.res.sum[, 'nfusion' := .N, by = bc]
  long.res.sum = long.res.sum[!duplicated(long.res.sum)]
  
  ## update detect infromation into seurat
  fusion_name = fusionParters[sampleName]
  
  seurat.rna$detected_fusion <- ifelse(colnames(seurat.rna) %in% long.res.sum$bc,
                                       fusion_name, seurat.rna$detected_fusion)
  seurat.rna$detected_nfusion[long.res.sum$bc] = long.res.sum$nfusion
  return(seurat.rna)
}


Add_wtInf2Seurat <- function(sampleName, seurat.rna, nmismatch = 0,
                                 dtype = 'both'){
  sampleName1 = sampleName
  if(sampleName %in% c('PAZBLA', 'PAYLNH', 'PAYUZJ')) 
    sampleName1 = paste0(sampleName, '_Redo')
  if(sampleName %in% c('HD2111_CD34', 'HD2689_Live'))
    sampleName1 = unlist(strsplit(sampleName, '_'))[1]
  
  message(paste0(sampleName, ':'))
  
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/LongRead/Final_Results/'
  read_file = paste0(dir0, sampleName1, '/', sampleName1, 
                     '.fastq.gz.demultiplexing.PASS.reads.txt')
  
  
  ## read all long reads for used bc
  long.res = fread(read_file)
  names(long.res)[c(1:3, 7)] = c('readname', 'bc', 'mismatch', 'sec_mismatch')
  long.res = long.res[mismatch <= nmismatch]
  long.res = subset(long.res, select = c('readname', 'bc', 'mismatch', 'sec_mismatch'))
  #long.res = long.res[sec_mismatch > 4]
  
  n0 = nrow(long.res)
  ## filtered out reads that are not map2 KMT2A PCR region
  if(sampleName == 'PAYKGI') return(seurat.rna)
  reads.pcr = fread(paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_reads2PCR.txt'))
  long.res = long.res[readname %in% reads.pcr$readname, ]
  #message(paste0(round(nrow(long.res)/n0 * 100, 1), '% reads in PCR!'))
  
  
  long.res$bc = paste0(sampleName, '_', long.res$bc)
  
  long.res = long.res[bc %in% colnames(seurat.rna)]
  
  
  ## read wt reads from alignqc
  if(dtype != 'longGF'){
    wt.file.alignqc = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_wt.txt')
    if(!file.exists(wt.file.alignqc)) return(seurat.rna)
    wt_reads.alignqc =  unique(fread(wt.file.alignqc)$readname)
    long.wt.alignqc = long.res[readname %in% wt_reads.alignqc & bc %in% colnames(seurat.rna)]
    
    message(paste('AlignQC detected', nrow(long.wt.alignqc), 'wt reads in',
                  length(unique(long.wt.alignqc$bc)), 'cells!'))
    
  }
  
  
  ## read wt reads from longGF
  if(dtype != 'alignqc'){
    wt_reads.longgf = fread(paste0(dir0, sampleName1, '/', 
                                   sampleName1, '.fastq.wt_reads.list'), header = F)$V1
    
    long.wt.longgf = long.res[readname %in% wt_reads.longgf & bc %in% colnames(seurat.rna)]
    
    message(paste('longGF detected', nrow(long.wt.longgf), 'wt reads in',
                  length(unique(long.wt.longgf$bc)), 'cells!'))
    
  }
  
  if(dtype == 'alignqc') long.wt = long.wt.alignqc
  if(dtype == 'longGF')  long.wt = long.wt.longgf
  
  if(dtype == 'both') {
    long.wt = long.wt.alignqc[readname %in% long.wt.longgf$readname]
    message(paste('Detected', nrow(long.wt), 'wt reads in',
                  length(unique(long.wt$bc)), 'cells by both methods!'))
    
  }
  
  message('-----------------------------------')
  if(nrow(long.wt) == 0) return(seurat.rna)
  
  
  long.res.sum = subset(long.wt, select = c(bc))
  long.res.sum[, 'nwt' := .N, by = bc]
  long.res.sum = long.res.sum[!duplicated(long.res.sum)]
  
  ## update detect infromation into seurat
  
  seurat.rna$detected_wt <- ifelse(colnames(seurat.rna) %in% long.res.sum$bc,
                                   'WT', seurat.rna$detected_wt)
  seurat.rna$detected_nwt[long.res.sum$bc] = long.res.sum$nwt
  
  
  ## write final read inf
  mdata0 = data.table(seurat.rna@meta.data, keep.rownames = T)
  mdata0$Ctype0[mdata0$Ctype0 %in% c('T', 'NKT')] = 'T/NK'
  setkey(mdata0, rn)
  long.wt[, 'Ctype0' := mdata0[J(long.wt$bc)]$Ctype0]
  write.table(subset(long.wt, select = c('readname', 'bc', 'Ctype0')),
              file = paste0('MetaData/scRNA/LongRead/Final_Fusions/', sampleName, '_wt.txt'),
              row.names = F, quote = F, sep = '\t')
  
  return(seurat.rna)
}

fusionRate_HD <- function(sampleName,  dtype = 'both'){
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/LongRead/Final_Results/'
  sampleName1 = unlist(strsplit(sampleName, '_'))[1]
  read_file = paste0(dir0, sampleName1, '/', sampleName1, 
                     '.fastq.gz.demultiplexing.PASS.reads.txt')
  
  
  ## read all long reads for used bc
  long.res = fread(read_file)
  names(long.res)[c(1:3, 7)] = c('readname', 'bc', 'mismatch', 'sec_mismatch')
  long.res = long.res[mismatch <= nmismatch]
  long.res = subset(long.res, select = c('readname', 'bc', 'mismatch', 'sec_mismatch'))
  #long.res = long.res[sec_mismatch > 4]
  
  ## filtered out reads that are not map2 KMT2A PCR region
  reads.pcr = fread(paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_reads2PCR.txt'))
  long.res = long.res[readname %in% reads.pcr$readname, ]
  message(paste0(round(nrow(long.res)/n0 * 100, 1), '% reads in PCR!'))
  
  
  message(paste0(sampleName, ':'))
  
  ## read fusion reads from longGF
  if(dtype != 'alignqc'){
    fusion.files = dir(paste0(dir0, sampleName1))
    fusion.file.sele = fusion.files[grepl(fusion.files, pattern = 'KMT2A|AP001267')]
    fusion.file.sele = fusion.file.sele[grepl(fusion.file.sele, pattern = 'list$')]
    fusion.file.sele = fusion.file.sele[grepl(fusion.file.sele, pattern = 'readname', fixed = T)]
    fusion_reads.longgf = lapply(fusion.file.sele, function(x) fread(paste0(dir0, sampleName1, '/', 
                                                                            x), header = F)$V1)
    fusion_reads.longgf = unique(do.call('c', fusion_reads.longgf))
    long.fusion.longgf = long.res[readname %in% fusion_reads.longgf ]
    
    message(paste('longGF detected', nrow(long.fusion.longgf), 'fusion reads in',
                  length(unique(long.fusion.longgf$bc)), 'cells!'))
    
  }
  
  ## read fusion reads from alignqc
  if(dtype != 'longGF'){
    fusion.file.alignqc = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
    if(!file.exists(fusion.file.alignqc)) return(seurat.rna)
    fusion_reads.alignqc =  unique(fread(fusion.file.alignqc)$readname)
    long.fusion.alignqc = long.res[readname %in% fusion_reads.alignqc]
    
    message(paste('AlignQC detected', nrow(long.fusion.alignqc), 'fusion reads in',
                  length(unique(long.fusion.alignqc$bc)), 'cells!'))
    
  }
  
  if(dtype == 'alignqc') long.fusion = long.fusion.alignqc
  if(dtype == 'longGF')  long.fusion = long.fusion.longgf
  
  if(dtype == 'both') {
    long.fusion = long.fusion.alignqc[readname %in% long.fusion.longgf$readname]
    message(paste('Detected', nrow(long.fusion), 'fusion reads in',
                  length(unique(long.fusion$bc)), 'cells by both methods!'))
    
  }
  return(unique(long.fusion$bc))
}

wtRate_HD <- function(sampleName,  dtype = 'both'){
  dir0 = '/mnt/isilon/tan_lab/chenc6/MLLr_Project/scRNA/LongRead/Final_Results/'
  sampleName1 = unlist(strsplit(sampleName, '_'))[1]
  read_file = paste0(dir0, sampleName1, '/', sampleName1, 
                     '.fastq.gz.demultiplexing.PASS.reads.txt')
  
  
  ## read all long reads for used bc
  long.res = fread(read_file)
  names(long.res)[c(1:3, 7)] = c('readname', 'bc', 'mismatch', 'sec_mismatch')
  long.res = long.res[mismatch <= nmismatch]
  long.res = subset(long.res, select = c('readname', 'bc', 'mismatch', 'sec_mismatch'))
  #long.res = long.res[sec_mismatch > 4]
  
  ## filtered out reads that are not map2 KMT2A PCR region
  if(sampleName == 'PAYKGI') return(seurat.rna)
  reads.pcr = fread(paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_reads2PCR.txt'))
  long.res = long.res[readname %in% reads.pcr$readname, ]
  message(paste0(round(nrow(long.res)/n0 * 100, 1), '% reads in PCR!'))
  
  message(paste0(sampleName, ':'))
  
  if(dtype != 'longGF'){
    wt.file.alignqc = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_wt.txt')
    if(!file.exists(wt.file.alignqc)) return(seurat.rna)
    wt_reads.alignqc =  unique(fread(wt.file.alignqc)$readname)
    long.wt.alignqc = long.res[readname %in% wt_reads.alignqc ]
    
    message(paste('AlignQC detected', nrow(long.wt.alignqc), 'wt reads in',
                  length(unique(long.wt.alignqc$bc)), 'cells!'))
    
  }
  
  
  ## read wt reads from longGF
  if(dtype != 'alignqc'){
    wt_reads.longgf = fread(paste0(dir0, sampleName1, '/', 
                                   sampleName1, '.fastq.wt_reads.list'), header = F)$V1
    
    long.wt.longgf = long.res[readname %in% wt_reads.longgf]
    
    message(paste('longGF detected', nrow(long.wt.longgf), 'wt reads in',
                  length(unique(long.wt.longgf$bc)), 'cells!'))
    
  }
  
  if(dtype == 'alignqc') long.wt = long.wt.alignqc
  if(dtype == 'longGF')  long.wt = long.wt.longgf
  
  if(dtype == 'both') {
    long.wt = long.wt.alignqc[readname %in% long.wt.longgf$readname]
    message(paste('Detected', nrow(long.wt), 'wt reads in',
                  length(unique(long.wt$bc)), 'cells by both methods!'))
    
  }
  
  return(unique(long.wt$bc))
}


plot_detectionRate <- function(mdata, do.test = F, test.bg = 'Mature_B',
                               type.read = 'fusion'){
  
  mdata$Ctype0[mdata$Ctype0 %in% c('T', 'NKT')] = 'T/NK'
  n_ctype = table(mdata$Ctype0)
  
  mdata0 = mdata[mdata$detected_nfusion > 0, ]
  if(type.read == 'wt') mdata0 = mdata[mdata$detected_nwt > 0, ]
  
  if(is.null(mdata0) | nrow(mdata0) < 1) return(NULL)
  n_detected = table(mdata0$Ctype0)
  
  ctype.not.detected = setdiff(names(n_ctype), names(n_detected))
  n_detected[ctype.not.detected] = 0
  
  n_ctype = n_ctype[names(n_detected)]
  
  
  data_fusion_ctype = data.table('Ctype' = names(n_ctype),
                                 'nCtype' = as.vector(n_ctype),
                                 'nDetected' = as.vector(n_detected))
  data_fusion_ctype[, 'rate' := nDetected/nCtype]
  
  data_fusion_ctype = data_fusion_ctype[Ctype != 'pDC']
  
  pvs = rep(1, nrow(data_fusion_ctype))
  names(pvs) = data_fusion_ctype$Ctype
  if(do.test){
    for(type0 in data_fusion_ctype$Ctype){
      if(type0 == test.bg) next
      pvs[type0] = pbinom(data_fusion_ctype[Ctype == type0]$nDetected, 
                          data_fusion_ctype[Ctype == type0]$nCtype,
                          prob = data_fusion_ctype[Ctype == test.bg]$rate,
                          lower.tail = F)
      pvs[type0] = min(pvs[type0], 1-pvs[type0])
    }
    data_fusion_ctype$pv = paste0('p=', format(pvs, digit = 2))
    data_fusion_ctype[Ctype == test.bg]$pv = 'background' 
  }
  
  data_fusion_ctype$rate = round(data_fusion_ctype$rate, 4)
  
  data_fusion_ctype0 = data_fusion_ctype
  
  #data_fusion_ctype$Ctype = factor(data_fusion_ctype$Ctype,
  #                                 levels = c('Blasts', 'Progenitors',
  #                                            'T/NK', 'Monocytes', 'Mature_B'))
  p1 <- ggplot(data = data_fusion_ctype, aes(x = Ctype, y = rate, fill = Ctype)) +
    geom_bar(stat = 'identity') + NoLegend() +
    xlab('') + ylab('Fraction') 
  
  if(do.test) p1 = p1 + geom_text(aes(label = pv, x = Ctype, y = rate+0.001))
   
  return(list('plot' = p1, 'pdata' = data_fusion_ctype0))
}

## add fusion inf to seurat  ####

seurat.rna = readRDS(file = 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds' )


sampleNames = c("PAYYBG", "PAZBSZ", "PAZFPH", "PAZGKI",
                "PAYWKL", "PAYYNY", "PAYZWN ", "PAYUZM",
                "154",   "PAYSBA", "PAYWJZ", "PAZBLA",
                "PAYKGI", "PAYLNH", "PAYUZJ", "PAYZLC",
                "PAYZVY", "PAZBGV")

fusionParters = c('MLLT1', 'MLLT1', 'MLLT1', 'MLLT3', 
                  'AFF1', 'AFF1', 'AFF1', 'MLLT10', 
                  'MLLT3', 'AFF1', 'MLLT1', 'MLLT10',
                  'AFF1', 'EPS15', 'MLLT1', 'EPS15',
                  'PSMF1', 'MLLT3')
names(fusionParters) = sampleNames

seurat.rna$detected_wt = 'NA'
seurat.rna$detected_nwt = 0
seurat.rna$detected_fusion = 'NA'
seurat.rna$detected_nfusion = 0
seurat.rna$nLongread = 0

detect_method = 'both'
for(name0 in sampleNames){
  seurat.rna <- Add_fusionInf2Seurat(name0, seurat.rna, 0, detect_method)
  seurat.rna <- Add_wtInf2Seurat(name0, seurat.rna, 0, detect_method)
}

#saveRDS(seurat.rna, file = 'Seurat_Objects/scRNA/seurat_regrCycleHeatShockGenes_pool_18Infants_scRNA_VEG3000_updated.rds' )
mdata.all = data.table(seurat.rna@meta.data, keep.rownames = T)
mdata.all = mdata.all[Ctype_Final %notin% c('doublets', 'Unknown')]
saveRDS(mdata.all, file='MetaData/scRNA/LongRead/mdata_with_longread_inf.rds')


## for HD ####
seurat.hd <- readRDS('Seurat_Objects/scRNA/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds')
seurat.hd = subset(seurat.hd, sample %in% c('donor2111_CD34p', 'donor2689_live'))
cnames = colnames(seurat.hd)
cnames = sapply(cnames, function(x) gsub('donor2111_CD34p', 'HD2111_CD34', x))
cnames = sapply(cnames, function(x) gsub('donor2689_live', 'HD2689_Live', x))
names(cnames) = NULL

samples = ifelse(grepl(cnames, pattern = 'HD2111_CD34'), 'HD2111_CD34', 'HD2689_Live')
seurat.hd$sample = samples
seurat.hd = RenameCells(seurat.hd, new.names = cnames)

seurat.hd$detected_fusion = 'NA'
seurat.hd$detected_nfusion = 0
seurat.hd$detected_wt = 'NA'
seurat.hd$detected_nwt = 0
seurat.hd$nLongread = 0

fusionParters[c('HD2111_CD34', 'HD2689_Live')] = 'DETECTED'

seurat.hd <- Add_fusionInf2Seurat('HD2111_CD34', seurat.hd, 0, 'both')
seurat.hd <- Add_wtInf2Seurat('HD2111_CD34', seurat.hd, 0, 'both')

seurat.hd <- Add_fusionInf2Seurat('HD2689_Live', seurat.hd, 0, 'both')
seurat.hd <- Add_wtInf2Seurat('HD2689_Live', seurat.hd, 0, 'both')
saveRDS(seurat.hd, file = 'Seurat_Objects/scRNA/hd_wLongread.rds')

mdata.hd = data.table(seurat.hd@meta.data, keep.rownames = T)
mdata.hd[, 'nLongread_str' := detected_nfusion + detected_nwt]

mdata.hd.pcr = mdata.hd[nLongread_str >= 1]

f2111 <- nrow(mdata.hd.pcr[detected_nfusion > 0 & sample == 'HD2111_CD34'])
w2111 <- nrow(mdata.hd.pcr[detected_nwt > 0 & sample == 'HD2111_CD34'])
f2689 <- nrow(mdata.hd.pcr[detected_nfusion > 0 & sample == 'HD2689_Live'])
w2689 <- nrow(mdata.hd.pcr[detected_nwt > 0 & sample == 'HD2689_Live'])


n2111 = nrow(mdata.hd.pcr[sample == 'HD2111_CD34'])
n2689 = nrow(mdata.hd.pcr[sample == 'HD2689_Live'])

c(f2111, w2111)/n2111
c(f2689, w2689)/n2689

hdata <- data.table('perc' = c(c(f2111, w2111)/n2111, c(f2689, w2689)/n2689) *100,
                    'sample' = rep(c('2111_CD34', '2689_Live'), each=2),
                    'read.type' = rep(c('Fusion', 'WT'), 2))
hdata$perc = round(hdata$perc, 2)
ggplot(hdata, aes(x=sample, y=perc, fill = read.type)) +
  geom_bar(stat = 'identity', position='dodge') + xlab('') +
  ylab('Detection rate: Percentage') + 
  geom_text(aes(label=paste0(perc, '%')), vjust=-0.5,  
            size=3.5, hjust = 0.5) + theme_classic() +
  theme(legend.title = element_blank())




## pool all data ####
mdata.all = readRDS(file='MetaData/scRNA/LongRead/mdata_with_longread_inf.rds')
mdata.all$Ctype0[mdata.all$Ctype0 %in% c('T', 'NKT')] = 'T/NK'
mdata.all[, 'nLongread_str' := detected_nfusion + detected_nwt]
mdata.all = mdata.all[sample != 'PAZGKI']
mdata.working = mdata.all[nLongread_str > 0]
#mdata.working = mdata.all[nLongread > 0]
res1 = plot_detectionRate(mdata.working, do.test = F, type.read = 'fusion')
res2 = plot_detectionRate(mdata.working, do.test = F, type.read = 'wt')

dd1 = res1$pdata
dd2 = res2$pdata
dd1$read.type = 'Fusion'
dd2$read.type = 'WT'
dd_fusion_wt <- rbind(dd1, dd2)  
pdata = subset(dd_fusion_wt, select = c('Ctype', 'rate', 'read.type'))
data.comb = pdata
data.comb$perc = round(data.comb$rate, 4) * 100
data.comb = data.comb[Ctype != 'Normal']
#data.comb$Ctype = factor(data.comb$Ctype, levels = c('Blasts', 'Progenitors',
#                                                     '2111_CD34', '2689_Live'))
p0 <- ggplot(data = data.comb, aes(x = Ctype, y = rate, fill = read.type)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  scale_fill_manual(values = brewer.pal(5, name = 'Set1')[1:2]) +
  geom_text(aes(label=paste0(perc, '%')), vjust=-0.5,  
            size=3.5, hjust = 0.5) +
  theme_classic() 


## combine pt and hd data
pdata = subset(dd_fusion_wt, select = c('Ctype', 'rate', 'read.type'))
hdata = subset(hdata, select = c('sample', 'perc', 'read.type'))
names(hdata) = names(pdata)
hdata$rate = hdata$rate/100
data.comb = rbind(pdata, hdata)
data.comb$perc = round(data.comb$rate, 4) * 100
data.comb = data.comb[Ctype != 'Normal']
data.comb$Ctype = factor(data.comb$Ctype, levels = c('Blasts', 'Progenitors',
                                                     '2111_CD34', '2689_Live'))
p7 <- ggplot(data = data.comb, aes(x = Ctype, y = rate, fill = read.type)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  scale_fill_manual(values = brewer.pal(5, name = 'Set1')[1:2]) +
  geom_text(aes(label=paste0(perc, '%')), vjust=-0.5,  
            size=3.5, hjust = 0.5) +
  theme_classic() 

ggsave(p7, filename = paste0('Figures/summary/LongRead/ExactMatch/fusionANDwt_rate_barplot_by_', 
                             detect_method, '.eps'), 
       device = 'eps', width = 6, height = 6)


## using all cells as denominator -- final ####
seurat.hd <- readRDS(file = 'Seurat_Objects/scRNA/hd_wLongread.rds')
mdata.hd = data.table(seurat.hd@meta.data, keep.rownames = T)
mdata.hd.pcr = mdata.hd[nLongread >= 0]

f2111 <- nrow(mdata.hd.pcr[detected_nfusion > 0 & sample == 'HD2111_CD34'])
w2111 <- nrow(mdata.hd.pcr[detected_nwt > 0 & sample == 'HD2111_CD34'])
f2689 <- nrow(mdata.hd.pcr[detected_nfusion > 0 & sample == 'HD2689_Live'])
w2689 <- nrow(mdata.hd.pcr[detected_nwt > 0 & sample == 'HD2689_Live'])


n2111 = nrow(mdata.hd.pcr[sample == 'HD2111_CD34'])
n2689 = nrow(mdata.hd.pcr[sample == 'HD2689_Live'])

c(f2111, w2111)/n2111
c(f2689, w2689)/n2689

hdata <- data.table('rate' = c(c(f2111, w2111)/n2111, c(f2689, w2689)/n2689) ,
                    'Ctype' = rep(c('2111_CD34', '2689_Live'), each=2),
                    'read.type' = rep(c('Fusion', 'WT'), 2))
hdata = subset(hdata, select = c(Ctype, rate, read.type))



mdata.all = readRDS(file='MetaData/scRNA/LongRead/mdata_with_longread_inf.rds')
mdata.all$Ctype0[mdata.all$Ctype0 %in% c('T', 'NKT')] = 'T/NK'
mdata.all[, 'nLongread_str' := detected_nfusion + detected_nwt]
mdata.all = mdata.all[sample != 'PAZGKI']

progenitor.g <- c('PAZBLA', 'PAYUZJ', 'PAYLNH',
                  '1154', 'PAYWKL', 'PAYKGI', 
                  'PAYUZM', 'PAYSBA', 'PAZGKI')
nprog.g = setdiff(sampleNames, progenitor.g)

mdata.working = mdata.all

res1 = plot_detectionRate(mdata.working, do.test = F, type.read = 'fusion')
res2 = plot_detectionRate(mdata.working, do.test = F, type.read = 'wt')

dd1 = res1$pdata
dd2 = res2$pdata
dd1$read.type = 'Fusion'
dd2$read.type = 'WT'
dd_fusion_wt <- rbind(dd1, dd2)  
pdata = subset(dd_fusion_wt, select = c('Ctype', 'rate', 'read.type'))

data.comb = rbind(pdata, hdata)
data.comb$perc = round(data.comb$rate, 4) * 100


data.comb$Ctype = factor(data.comb$Ctype, levels = c('Blasts', 'Progenitors',
                                                     'Mature_B', 'Monocytes', 'T/NK',
                                                     '2111_CD34', '2689_Live'))
p1 <- ggplot(data = data.comb, aes(x = Ctype, y = rate, fill = read.type)) +
  geom_bar(stat = 'identity', position='dodge') + 
  xlab('') + ylab('Fraction') + 
  scale_fill_manual(values = brewer.pal(5, name = 'Set1')[1:2]) +
  geom_text(aes(label=paste0(perc, '%')), vjust=-0.5,  
            size=3.5, hjust = 0.5) +
  theme_classic() 

ggsave(p1, filename = 'Figures/summary/LongRead/ExactMatch/updated_frac_fusion_wt_allCells.eps',
       device = 'eps', width = 6, height = 6)
