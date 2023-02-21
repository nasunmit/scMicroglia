
#####
##########################
id='iMGs.ctrl_LPS.snRNAseq.May2022.seurat'
##########################
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
wr <- colorRampPalette(colors = c( "white", "red"))(100)
rwb <- colorRampPalette(colors = c("blue", "white", "red"))(100)
ryb=colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)
prgn=colorRampPalette(rev(brewer.pal(n = 10, name ="PRGn")))(100)
wb <- colorRampPalette(colors = c( "white", "blue"))(100)
#############################
brain=readRDS('../iMGs.snRNAseq.May2022.seurat.rds')
table(brain$condition)
meta=brain@meta.data
meta=meta[meta$condition %in% c('Ctrl','LPS'),]
dim(meta)
brain=subset(brain,cells=rownames(meta))

brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
brain <- ScaleData(brain,features =VariableFeatures(brain), verbose=FALSE)
brain <- RunPCA(brain,features = VariableFeatures(brain), npcs = 50,verbose=FALSE)
ElbowPlot(brain,ndims = 50)
k=1:30
brain <- RunUMAP(brain, reduction = "pca", dims = k, verbose=FALSE)
DimPlot(brain)
brain <- FindNeighbors(brain, dims = k)
brain <- FindClusters(brain, resolution = 0.6)


pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T)
DimPlot(brain, reduction = "umap",label=F)
DimPlot(brain,group.by = 'library')
DimPlot(brain,group.by = 'sample')
DimPlot(brain,group.by = 'condition')
dev.off()

pdf(file=paste(id,'.umap_split.pdf',sep=''),width = 10,height = 5)
DimPlot(brain, reduction = "umap",label=T,split.by = 'condition')
DimPlot(brain, reduction = "umap",label=F,split.by = 'condition')
dev.off()

saveRDS(brain,file=paste(id,'rds',sep='.'))

## Finding differentially expressed genes (cluster biomarkers)
markers <- FindAllMarkers(object = brain, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(markers,file=paste(id,'.clusterDEGs.txt',sep = ''),sep='\t',quote=F)
table(markers$cluster)


##### 
id=paste(id,'.harmony',sep='')
brain <- RunHarmony(brain, group.by.vars = "library")
brain <- RunUMAP(brain, reduction = "harmony", dims = 1:30)
brain <- FindNeighbors(brain, reduction = "harmony", dims = 1:30) %>% FindClusters()
DimPlot(brain, group.by = c("sample", "ident", "condition"), ncol = 3)

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T)
DimPlot(brain, reduction = "umap",label=F)
DimPlot(brain,group.by = 'library')
DimPlot(brain,group.by = 'sample')
DimPlot(brain,group.by = 'condition')
dev.off()

pdf(file=paste(id,'.umap_split.pdf',sep=''),width = 10,height = 5)
DimPlot(brain, reduction = "umap",label=T,split.by = 'condition')
DimPlot(brain, reduction = "umap",label=F,split.by = 'condition')
dev.off()

## Finding differentially expressed genes (cluster biomarkers)
markers <- FindAllMarkers(object = brain, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(markers,file=paste(id,'.clusterDEGs.txt',sep = ''),sep='\t',quote=F)
table(markers$cluster)

saveRDS(brain,file=paste(id,'rds',sep='.'))



