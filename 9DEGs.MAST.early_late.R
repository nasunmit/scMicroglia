##########################
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
###########################

## Identify DEGs
#######
library(edgeR)
library(limma)
library(MAST)
library(tidyr)
######
f='ROSMAP.Microglia.6regions.seurat.harmony.selected.rds'
####

fid=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]rds')[[1]][1]
brain=readRDS(f)
meta=brain@meta.data
table(meta$ADdiag3types)
metasel=meta[meta$ADdiag3types!='nonAD',]
brain=subset(brain,cells=rownames(metasel))
brain$cellsubtype=paste('MG',brain$seurat_clusters,sep='')
meta=brain@meta.data
data=brain@assays$RNA@counts
lib_size <- colSums(data)
norm <- t(t(data)/lib_size * median(lib_size))
## cell type
tb=table(meta$cellsubtype)
tb=tb[tb>100]
ctps=names(sort(tb,decreasing = T))

for (c in ctps){
  tmpmeta=meta[meta$cellsubtype==c,]
  tmpmeta=droplevels(tmpmeta)
  tmpdata=data[,rownames(tmpmeta)]
  tmpnorm=norm[,rownames(tmpmeta)]
  group=as.character(tmpmeta$ADdiag3types)
  table(group)
  
  bregion=tmpmeta$brainRegion
  nfeat=tmpmeta$nFeature_RNA
  batch=tmpmeta$batch
  dxpark=tmpmeta$dxpark
  dlb=tmpmeta$clin_dlb
  age=tmpmeta$age_death
  sex=tmpmeta$msex
  pmi=tmpmeta$pmi
  mt=tmpmeta$percent.mt
  rp=tmpmeta$percent.rp
  #### MAST   age+sex+pmi+tdp
  log_counts <- log2(tmpdata + 1)
  fData <- data.frame(names = rownames(log_counts))
  rownames(fData) <- rownames(log_counts)
  cData <- data.frame(cond = group,br=bregion,batch=batch,dxpark=dxpark,dlb=dlb,age=age,sex=sex,pmi=pmi,mt=mt,rp=rp)
  cData$cond=relevel(factor(cData$cond),"earlyAD")
  rownames(cData) <- colnames(log_counts)
  obj <- FromMatrix(as.matrix(log_counts), cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
  cond <- factor(colData(obj)$cond)
  if (length(table(batch))>1){
    zlmCond <- zlm(~ cond + cngeneson+br+batch+dxpark+dlb+age+sex+pmi+mt+rp, obj)
  }else{zlmCond <- zlm(~ cond + cngeneson+br+dxpark+dlb+age+sex+pmi+mt+rp, obj) }
  summaryCond <- summary(zlmCond, doLRT = "condlateAD")
  
  summaryDt <- summaryCond$datatable
  res <- merge(summaryDt[contrast=='condlateAD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
               summaryDt[contrast=='condlateAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  
  res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  res=as.data.frame(res)
  
  write.table(res,file=paste(paste(fid,c,sep = '.'),'.early_late.MASTres.txt',sep=''),sep = '\t',row.names = F,quote = F)
}


