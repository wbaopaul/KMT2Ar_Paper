source('scDataAnalysis_Utilities.R')
`%notin%` = Negate(`%in%`) 
score.thr = 2

## for AF9 (MLLT3) ####
for(sampleName in c('MLLr1154', 'MLLr882304', 'MLLr879339')){
  dd = fread(paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz'), 
             skip = 1)
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
  
  if(sampleName != 'MLLr1154') dd[, 'chr' := paste0('chr', chr)]
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr9' & start >= 20341669 & end <= 20622499 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}


## for AF4 ####
for(sampleName in c('MLLr876545', 'MLLr877780', 'MLLr874013', 'MLLr870684', 'MLLr878516')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
  dd[, 'chr' := paste0('chr', chr)]
  score.thr = 2
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr4' & start >= 87006988 & end <= 87141038 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}

## for ENL ####
for(sampleName in c('MLLr877476', 'MLLr879583', 'MLLr881823', 'MLLr876533', 'MLLr875703')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
  if(sampleName %notin% c('MLLr875703')) dd[, 'chr' := paste0('chr', chr)]
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr19' & start >= 6210381 & end <= 6279975 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}

## for AF10 ####
for(sampleName in c( 'MLLr879440', 'MLLr875706')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
 
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr10' & start >= 21534232 & end <= 21743630 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}

## for EPS15 ####
for(sampleName in c( 'MLLr871427', 'MLLr878289')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
  if(sampleName != 'MLLr871427') dd[, 'chr' := paste0('chr', chr)]

  
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr1' & start >= 51354263 & end <= 51519266 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}

## for PSMF1 ####
for(sampleName in c( 'MLLr878501')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]
  dd[, 'chr' := paste0('chr', chr)]
  
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr == 'chr20' & start >= 1113263 & end <= 1167782 & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}

## fusion for HD ####
for(sampleName in c( 'HD2111_CD34', 'HD2689_Live')){
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(chimera.file)) next 
  dd = fread(chimera.file, skip = 1)
  names(dd)[c(1:4, 10)] = c('chr', 'start', 'end', 'readname', 'score')
  dd[, 'len' := end-start]

  
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118481810 & start < 118481830) | (end >= 118481810 & end < 118481830)]
  
  
  dd2 = dd[ chr != 'chr11' & score >= score.thr]
  sele.reads = intersect(dd1$readname, dd2$readname)
  
  dd.sele = dd[readname %in% sele.reads]
  message(paste0(sampleName, ': detected ', length(sele.reads), ' fusion reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_fusion.txt')
  write.table(dd.sele, file = outFile,
              row.names = F, quote = F, sep = '\t')
  
}


## for WT ####
sampleNames = c("MLLr877476", "MLLr879583", "MLLr881823", "MLLr882304",
                "MLLr876545", "MLLr877780", "MLLr878516", "MLLr875706",
                "MLLr1154",   "MLLr874013", "MLLr876533", "MLLr879440",
                "MLLr870684", "MLLr871427", "MLLr875703", "MLLr878289",
                "MLLr878501", "MLLr879339", "HD2111_CD34", 'HD2689_Live')

for(sampleName in sampleNames){
  long.read.file = paste0('AlignQC/Chimera_Output/', sampleName, '_best.sorted.bed.gz')
  chimera.file = paste0('AlignQC/Chimera_Output/', sampleName, '_chimera.bed.gz')
  if(!file.exists(long.read.file)) next 
  dd.trans = fread(chimera.file, skip = 1, select = c(1:4, 10))
  names(dd.trans) = c('chr', 'start', 'end', 'readname', 'score')
  
  dd = fread(long.read.file, skip = 1, select = c(1:4, 10))
  names(dd) = c('chr', 'start', 'end', 'readname', 'score')
  
  dd[, 'len' := end-start]
  dd.trans[, 'len' :=end-start]
  if(sampleName %notin% c('MLLr1154', 'MLLr879440', 'MLLr871427', 'HD2111_CD34',
                      'HD2689_Live', 'MLLr875703', 'MLLr875706')) {
    dd[, chr := paste0('chr', chr)]
    dd.trans[, chr := paste0('chr', chr)]
  }
  dd0 = dd[chr == 'chr11']
  dd0 = dd0[(start >= 118481713 & start <= 118482093) | (end >= 118481713 & end <= 118482093)]
  
  dd0.trans = dd.trans[chr == 'chr11']
  dd0.trans = dd0.trans[(start >= 118481713 & start <= 118482093) | (end >= 118481713 & end <= 118482093)]
  dd.pcr = rbind(dd0, dd0.trans[dd0.trans$readname %notin% dd0$readname, ])
  
  
  outFile0 = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_reads2PCR.txt')
  write.table(dd.pcr, file = outFile0,
              row.names = F, quote = F, sep = '\t')
  
  dd1 = dd[chr == 'chr11' & score >= score.thr]
  dd1 = dd1[(start >= 118490129 & start < 118523917) | (end >= 118490129 & end < 118523917)]
  
  dd2 = dd[chr != 'chr11' & score >= score.thr]
  dd1 = dd1[readname %notin% dd2$readname]
  
  message(paste0(sampleName, ': detected ', nrow(dd1), ' wt reads!'))
  outFile = paste0('MetaData/scRNA/LongRead/', sampleName, '_alignqc_filtered_wt.txt')
  write.table(dd1, file = outFile,
              row.names = F, quote = F, sep = '\t')
}

