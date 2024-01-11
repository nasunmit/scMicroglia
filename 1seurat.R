setwd('~/Dropbox (Personal)/work/projects/project.brain/project0.microglia/analysis/decont_counts/1.mic/1seurat.subset/')

source('~/Dropbox (MIT)/work/dataset/Rscripts/rcolors.R')
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
###########################
id='ROSMAP.Microglia.6regions.seurat'

### data
brain=readRDS('~/Dropbox (Personal)/work/projects/project.brain/project0.microglia/analysis/decont_counts/0.immune/ROSMAP.ImmuneCells.6regions.seurat.harmony.rds')
meta=brain@meta.data ## 13: T, 15: neuron, 14: Oligo
clusters=c(0:12)
meta=meta[meta$seurat_clusters %in% clusters,]
dim(meta)
#### seurat
brain=subset(brain,cells=rownames(meta))
brain


brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)
## scaling the data
brain <- ScaleData(brain, features = rownames(brain))
## Perform linear dimensional reduction
brain <- RunPCA(brain, features = VariableFeatures(object = brain))
DimPlot(brain, reduction = "pca")
ElbowPlot(brain,ndims = 50)
k=1:30

## cluster the cells
brain <- FindNeighbors(brain, dims = k)
brain <- FindClusters(brain, resolution = 0.5)
head(Idents(brain), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
brain <- RunUMAP(brain, dims = k)
mycol=col.all[1:length(table(Idents(brain)))]

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T,cols=mycol,raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'batch',raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'brainRegion',cols=c(col.set1,col.set3),raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))),raster=FALSE)+NoLegend()
dev.off()

saveRDS(brain,file=paste(id,'.rds',sep=''))
write.table(brain@meta.data,file=paste(id,'.meta.txt',sep=''),sep = '\t',quote = F)
############
##### batch correctio: run harmony
id=paste(id,'.harmony',sep='')
options(repr.plot.height = 2.5, repr.plot.width = 6)
#brain <- brain %>% RunHarmony("batch", plot_convergence = TRUE)
#brain <- brain %>% RunHarmony("library_id", plot_convergence = TRUE)
brain <- brain %>% RunHarmony("projid", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(brain, 'harmony')
harmony_embeddings[1:5, 1:5]

brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

mycol=col.all[1:length(table(Idents(brain)))]
DimPlot(brain, reduction = "umap",label=T,cols=mycol)
DimPlot(brain, reduction = "umap",label=F,group.by = 'batch')
DimPlot(brain, reduction = "umap",label=F,group.by = 'brainRegion',cols=c(col.set1,col.set3))
DimPlot(brain, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))))+NoLegend()


meta=brain@meta.data
pdf(file=paste(id,'.umap_eachcluster.pdf',sep=''),width = 6,height = 5)
for (i in names(table(meta$seurat_clusters))){
  print(DimPlot(brain, reduction = "umap",label=T,cells=rownames(meta[meta$seurat_clusters==i,])))
}
dev.off()

###

table(Idents(brain))
brain=subset(brain,idents=c(0:12))

####
require(spatstat)
require(sp)
require(raster)

umaps=brain@reductions$umap@cell.embeddings
cids=names(table(Idents(brain)))

pdf(file='distance_density.pdf')
for (i in cids){
  tmpmeta=meta[meta$seurat_clusters==i,]
  tmpumaps=umaps[rownames(tmpmeta),]
  x=tmpumaps[,1]
  y=tmpumaps[,2]
  p <- SpatialPoints(coords = matrix(c(x, y), ncol = 2))
  plot(p,main=i)
  pp = ppp(x,y, c(min(x),max(x)), c(min(y),max(y))) # all points in a (0,1) default window
  d <- density.ppp(pp, sigma = 0.1)
  dp <- density.ppp(pp, sigma = 0.1, at="points")
  dr = raster(d)
  mycenter=xyFromCell(dr, which.max(dr))
  mydist=pointDistance(mycenter,tmpumaps,lonlat = F)
  names(mydist)=rownames(tmpumaps)
 # mycells=names(mydist[mydist<4])
  print(plot(density(mydist),main=i))
  #plot(tmpumaps[mycells,1],tmpumaps[mycells,2])
  
}
dev.off()

cutoffs=c(4,5,4,5,4,4,3,4,3,2,3,4,2)
names(cutoffs)=cids
selcells=c()
for (i in cids){
  tmpmeta=meta[meta$seurat_clusters==i,]
  tmpumaps=umaps[rownames(tmpmeta),]
  x=tmpumaps[,1]
  y=tmpumaps[,2]
  p <- SpatialPoints(coords = matrix(c(x, y), ncol = 2))
  #plot(p,main=i)
  pp = ppp(x,y, c(min(x),max(x)), c(min(y),max(y))) # all points in a (0,1) default window
  d <- density.ppp(pp, sigma = 0.1)
  dp <- density.ppp(pp, sigma = 0.1, at="points")
  dr = raster(d)
  mycenter=xyFromCell(dr, which.max(dr))
  mydist=pointDistance(mycenter,tmpumaps,lonlat = F)
  names(mydist)=rownames(tmpumaps)
  mycells=names(mydist[mydist<cutoffs[i]])
  selcells=c(selcells,mycells)
  #print(plot(density(mydist),main=i))
  #plot(tmpumaps[mycells,1],tmpumaps[mycells,2])
  
}
length(selcells)

#############
id=paste(id,'selected',sep='.')
brainsel=subset(brain,cells=selcells)
brainsel=subset(brainsel,idents=c(0:8,10:12))
table(Idents(brainsel))
meta=brainsel@meta.data
meta=droplevels(meta)
brainsel@meta.data=meta

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brainsel, reduction = "umap",label=T,cols=mycol,raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'batch',raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'brainRegion',cols=c(col.set1,col.set3),raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))),raster=FALSE)+NoLegend()
dev.off()

saveRDS(brainsel,file=paste(id,'.rds',sep=''))
write.table(brainsel@meta.data,file=paste(id,'.meta.txt',sep=''),sep = '\t',quote = F)


# find markers for every cluster compared to all remaining cells, report only the positive ones
brain.markers <- FindAllMarkers(brainsel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(brain.markers$cluster)
top5=brain.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.table(brain.markers,file=paste(id,'.clusterDEGs.txt',sep=''),sep = '\t',quote = F)


###########
## pseudo-bulk on DEGs
newmeta=brainsel@meta.data
data=brainsel@assays$RNA@data[unique(brain.markers$gene),]
bulkdata=c()
myclusters=names(table(brain.markers$cluster))
for (c in myclusters){
  cells=rownames(newmeta[newmeta$seurat_clusters==c,])
  tmp=data[,cells]
  bulkdata=cbind(bulkdata,rowMeans(as.matrix(tmp)))
}
colnames(bulkdata)=paste('C',myclusters,sep='')
write.table(bulkdata,file=paste(id,'.clusterDEGs_pseudobulk.txt',sep=''),sep = '\t',quote = F)

##
## pseudo-bulk on all genes
newmeta=brainsel@meta.data
data=brainsel@assays$RNA@data
bulkdata=c()
myclusters=names(table(newmeta$seurat_clusters))
for (c in myclusters){
  cells=rownames(newmeta[newmeta$seurat_clusters==c,])
  tmp=data[,cells]
  bulkdata=cbind(bulkdata,rowMeans(as.matrix(tmp)))
}
colnames(bulkdata)=paste('MG',myclusters,sep='')
write.table(bulkdata,file=paste(id,'.cluster_pseudobulk.txt',sep=''),sep = '\t',quote = F)
###

### check known markers of cell types in brain 
knowngenes=as.character(read.table('~/Dropbox (MIT)/work/dataset/markers/human.brain.knownMarkers.txt',header = F,sep = '\t')$V1)
genes=intersect(rownames(brainsel),knowngenes)
pdf(file=paste(id,'.umap.brain_knownMarkers.pdf',sep=''),width = 6,height = 5)
for (g in genes){
  print(FeaturePlot(brainsel, features = g,raster=FALSE))
}
dev.off()

pdf(file=paste(id,'.dotplot.brain_knownMarkers.pdf',sep=''),width = 10,height = 7)
DotPlot(brainsel,features = genes)+RotatedAxis()
dev.off()

### check known markers of cell types in immune cells 
ref1=as.character(read.table('~/Dropbox (MIT)/work/dataset/microglia/microglia.markers.human.txt',header = F,sep = '\t')$V1)

ref2=read.table('~/Dropbox (MIT)/work/dataset/immune.signature/immune-cell-markers.list.txt',header = T,sep = '\t')
ref2=strsplit(paste(ref2$symbol,collapse = ','),',')[[1]]

ref3=as.character(read.table('~/Dropbox (MIT)/work/dataset/immune.signature/Human.PBMC_CSF_Microglia.panelGenes.txt',header = F,sep = '\t')$V1)

genes=unique(c(ref1,ref2,ref3))
genes=intersect(genes,rownames(brainsel))

pdf(file=paste(id,'.umap.immune_knownMarkers.pdf',sep=''),width = 6,height = 5)
for (g in genes){
  print(FeaturePlot(brainsel, features = g,raster=FALSE))
}
dev.off()

pdf(file=paste(id,'.dotplot.immune_knownMarkers.pdf',sep=''),width = 30,height = 7)
DotPlot(brainsel,features = genes)+RotatedAxis()
dev.off()


#### cell type score

ref=read.table('~/Dropbox (MIT)/work/dataset/immune.signature/human.22_immune_celltypes.signature.matrix.txt',header = T,sep = '\t')

###### cell type signature
n=ncol(brainsel@meta.data)
data=ref[rownames(ref) %in% brain.markers$gene,]
for (i in c(1:ncol(data))){
  genes=rownames(data[data[,i]==1,])
  myfeature=list(genes)
  if(length(genes)>1){
    brainsel=AddModuleScore(brainsel,features = myfeature,name=paste(colnames(data)[i],'.',sep=''))
  }
}

## rename meta for signature
head(brainsel@meta.data)
colnames(brainsel@meta.data)=c(colnames(brainsel@meta.data)[1:n],gsub('[.]1','',colnames(brainsel@meta.data)[(n+1):ncol(brainsel@meta.data)]))

pdf(file=paste(id,'.umap.immune_celltype_score.pdf',sep=''),width = 6,height = 5)
for (g in colnames(brainsel@meta.data)[(n+1):ncol(brainsel@meta.data)]){
  print(FeaturePlot(brainsel, features = g,cols = c('white','black'),raster=FALSE))
}
dev.off()

write.table(brainsel@meta.data,file=paste(id,'.metadata.immune_celltype_score.txt',sep=''),sep = '\t',quote=F)

myft=colnames(brainsel@meta.data)[(n+1):ncol(brainsel@meta.data)]
pdf(file=paste(id,'.dotplot.immune_celltype_score.pdf',sep=''),width = 20,height = 10)
print(DotPlot(brainsel,features = rev(myft))+ RotatedAxis())
dev.off()

####
pdf(file = paste(id,'.CAMs_signature.pdf',sep=''),width = 16,height = 14)
FeaturePlot(brainsel, features = c("MRC1", "LYVE1", "CD163", "SIGLEC1"),raster=FALSE) ## CAMs signature
dev.off()

pdf(file = paste(id,'.LymphocyteSignature.pdf',sep=''),width = 8,height = 7)
FeaturePlot(brainsel, features = c("TRAC", "TRBC2", "CD52", "IL32"),raster=FALSE)## Lymphocyte signature
dev.off()
pdf(file = paste(id,'.MyeloidSignature.pdf',sep=''),width = 8,height = 7)
FeaturePlot(brainsel, features = c("ITGAM", "MS4A6A", "TYROBP", "CD14"),raster=FALSE) ## Myeloid signature
dev.off()
pdf(file = paste(id,'.MonocyteSignature.pdf',sep=''),width = 8,height = 7)
FeaturePlot(brainsel, features = c("CCR2", "PLAC8", "CLEC12A", "FCN1"),raster=FALSE) ## Monocyte signature
dev.off()

####




