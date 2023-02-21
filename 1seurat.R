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
id='ROSMAP.ImmuneCells.6regions.seurat'

### count matrix and meta data 
data=readRDS('counts.rds')
meta=readRDS('meta.rds')
#### seurat
brain=CreateSeuratObject(counts = data, project = "immune",meta.data = meta)
brain
brain[["percent.rp"]] <- PercentageFeatureSet(brain, pattern = "^RP")
VlnPlot(brain, features = c("nFeature_RNA", "percent.rp", "percent.mt"), ncol = 3)

brain <- subset(brain, subset = percent.rp<5)
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
############
##### batch correctio: run harmony
id=paste(id,'.harmony',sep='')
options(repr.plot.height = 2.5, repr.plot.width = 6)
brain <- brain %>% RunHarmony("batch", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(brain, 'harmony')
harmony_embeddings[1:5, 1:5]

brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

table(Idents(brain))

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T,raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'batch',raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'brainRegion',raster=FALSE)
DimPlot(brain, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))),raster=FALSE)+NoLegend()
dev.off()

saveRDS(brain,file=paste(id,'.rds',sep=''))
write.table(brain@meta.data,file=paste(id,'.meta.txt',sep=''),sep = '\t',quote = F)


# find markers for every cluster compared to all remaining cells, report only the positive ones
brain.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(brain.markers$cluster)
write.table(brain.markers,file=paste(id,'.clusterDEGs.txt',sep=''),sep = '\t',quote = F)

knowngenes=c("SYT1","SNAP25","GRIN1","NRGN","SLC17A7","CAMK2A","GAD1","GAD2","AQP4","GFAP","MBP","MOBP","PLP1","CSF1R","CD74","C3","VCAN","PDGFRA","CSPG4","FLT1","CLDN5","AMBP")
othergenes=c('MRC1','CD163','LYVE1','ITGA4','KLPD1','CD3E','CD3G','CD44')
mygenes=intersect(c(knowngenes,othergenes),rownames(brain))
pdf(file='ROSMAP.ImmuneCells.6regions.seurat.harmony.dotplot.keygenes.pdf',width = 9,height = 4.5)
DotPlot(brain,features = mygenes,cluster.idents=T)
dev.off()

###########
## select microglia clusters ###
id='ROSMAP.Microglia.6regions.seurat'
brainsel=subset(brain,cells=selcells)
brainsel=subset(brainsel,idents=c(0:8,10:12)) #### 9: CAMs, 13: T, 15: neuron, 14: Oligo
table(Idents(brainsel))
meta=brainsel@meta.data
meta=droplevels(meta)
brainsel@meta.data=meta

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brainsel, reduction = "umap",label=T,raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'batch',raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'brainRegion',raster=FALSE)
DimPlot(brainsel, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))),raster=FALSE)+NoLegend()
dev.off()

saveRDS(brainsel,file=paste(id,'.rds',sep=''))
write.table(brainsel@meta.data,file=paste(id,'.meta.txt',sep=''),sep = '\t',quote = F)


# find markers for every cluster compared to all remaining cells, report only the positive ones
brain.markers <- FindAllMarkers(brainsel, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(brain.markers$cluster)
write.table(brain.markers,file=paste(id,'.clusterDEGs.txt',sep=''),sep = '\t',quote = F)


