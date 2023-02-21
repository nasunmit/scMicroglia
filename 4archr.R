library(ArchR)
########
df=read.table('/home/nasun/data-lab5/projects/project0.microglia/analysis.snATAC/All.ATAC.samp.info.txt',header = T,sep = '\t')
table(df$region)
df=df[df$region!='MB',]


inputFiles=list.files(path = '/home/nasun/data-lab5/projects/project0.microglia/old/data/ShareWithNa.Mic/fragments',pattern = 'tsv.gz',full.names = T)
sNames=c()
for (f in inputFiles){nm=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]')[[1]][1];sNames=c(sNames,nm)}
names(inputFiles)=sNames

ov=intersect(sNames,df$SampID)
df=df[df$SampID %in% ov,]
sNames=ov
inputFiles=inputFiles[ov]
## setup 
addArchRThreads(threads = 16)
addArchRGenome("hg38") 
####
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,sampleNames = sNames,filterTSS = 6, filterFrags = 1000, addTileMat = TRUE,addGeneScoreMat = TRUE)

### 
ArrowFiles=list.files(pattern = 'arrow')
ArrowFiles

brain1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "brain.microglia",
  copyArrows = F #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
brain1

paste0("Memory Size = ", round(object.size(brain1) / 10^6, 3), " MB")
getAvailableMatrices(brain1)

head(brain1$cellNames)
head(brain1$Sample)

rownames(df)=df$SampID
libs=as.matrix(unname(as.data.frame(strsplit(brain1$cellNames,'#'))[1,]))[1,]
mymeta=df[libs,]
#mymeta=data.frame(lapply(mymeta, as.character), stringsAsFactors=FALSE)
rownames(mymeta)=brain1$cellNames
for (i in c(2:ncol(mymeta))){
  brain1=addCellColData(ArchRProj = brain1,data=mymeta[,i],name = colnames(mymeta)[i],cells=brain1$cellNames,force = T)
}
str(brain1@cellColData)
#Plotting Sample Statistics from an ArchRProject
# Make a ridge plot for each sample for the TSS enrichment scores.
p1 <- plotGroups(
  ArchRProj = brain1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)


# Make a violin plot for each sample for the TSS enrichment scores.

p2 <- plotGroups(
  ArchRProj = brain1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)


#Make a ridge plot for each sample for the log10(unique nuclear fragments).
p3 <- plotGroups(
  ArchRProj = brain1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)



#Make a violin plot for each sample for the log10(unique nuclear fragments).

p4 <- plotGroups(
  ArchRProj = brain1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)


plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = brain1, addDOC = FALSE, width = 30, height = 30)


## Plotting Sample Fragment Size Distribution and TSS Enrichment Profiles.
p1 <- plotFragmentSizes(ArchRProj = brain1)
p2 <- plotTSSEnrichment(ArchRProj = brain1)

plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = brain1, addDOC = FALSE, width = 30, height = 30)

#################################################################
##### Dimensionality Reduction with ArchR########################
#################################################################
brain2=brain1
brain2 <- addIterativeLSI(
  ArchRProj = brain2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 10, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.5), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
)

## batch correction with Harmony
brain2 <- addHarmony(
  ArchRProj = brain2,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

#################################################################
##### Clustering ################################################
#################################################################
## 
brain2 <- addClusters(
  input = brain2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5
)

table(brain2$Clusters)
cM <- confusionMatrix(paste0(brain2$Clusters), paste0(brain2$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p, name = "samples_cluster.pheatmap.pdf", ArchRProj = brain2, addDOC = FALSE, width = 10, height = 8)

#################################################################
##### UMAP ######################################################
#################################################################
### UMAP

brain2 <- addUMAP(
  ArchRProj = brain2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "region", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "projid", embedding = "UMAP")


plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brain2, addDOC = FALSE, width = 6, height = 6)


## Dimensionality Reduction After Harmony

brain2 <- addUMAP(
  ArchRProj = brain2, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

p1 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = brain2, colorBy = "cellColData", name = "region", embedding = "UMAPHarmony")

plotPDF(p1,p2,p3, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = brain2, addDOC = FALSE, width = 6, height = 6)

saveRDS(brain2,file='brain.microglia/brain2.rds')
#########
brain=readRDS('brain.microglia/brain2.rds')
meta=brain@cellColData
table(meta$Clusters)

remove<-c("C1","C6") ## remove potential doublet clusters
keep<-!(brain$Clusters %in% remove)
sum(keep=="TRUE") ## the number of cell to be kept
brain=brain[keep,]

brain@projectMetadata$outputDirectory="brain.microglia.filter"
#################################################################
##### Dimensionality Reduction with ArchR########################
#################################################################

brain <- addIterativeLSI(
  ArchRProj = brain,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 10, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(1), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force=TRUE
)

## batch correction with Harmony
brain <- addHarmony(
  ArchRProj = brain,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

#################################################################
##### Clustering ################################################
#################################################################
## 
brain <- addClusters(
  input = brain,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1,
  force = TRUE
)

table(brain$Clusters)
cM <- confusionMatrix(paste0(brain$Clusters), paste0(brain$Sample))
cM

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(p, name = "samples_cluster.pheatmap.pdf", ArchRProj = brain, addDOC = FALSE, width = 10, height = 8)

#################################################################
##### UMAP ######################################################
#################################################################
### UMAP

brain <- addUMAP(
  ArchRProj = brain, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "region", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "projid", embedding = "UMAP")


plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6)


## Dimensionality Reduction After Harmony

brain <- addUMAP(
  ArchRProj = brain, 
  reducedDims = "Harmony", 
  name = "UMAPHarmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p2 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p3 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "region", embedding = "UMAPHarmony")

plotPDF(p1,p2,p3, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6)


brain2=brain
#################################################################
##### identify marker genes ############################################
##############################################################
### 
markersGS <- getMarkerFeatures(
  ArchRProj = brain2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList,file='brain.microglia.filter/markerList.rds')

df=read.table('microglia.markers.human.txt',sep = '\t')
markerGenes.mic = as.character(df$V1)
rnadeg=read.table('ROSMAP.Microglia.6regions.seurat.harmony.clusterDEGs.txt',header = T,sep = '\t')
markerGenes=unique(as.character(rnadeg$gene))

######
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)

##
heatmapGS.mic <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
  labelMarkers = markerGenes.mic,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)
############
heatmapGS.both <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = c(markerGenes,markerGenes.mic),
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS.both, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker-Heatmap", width = 20, height = 6, ArchRProj = brain2, addDOC = FALSE)

## overlap between snRNA subtype markers and snATAC subtype markers
atacgenes=colnames(heatmapGS@matrix)
ovgenes=intersect(atacgenes,markerGenes)
rnadeg.sel=rnadeg[rnadeg$gene %in% ovgenes,]

######### ovgenes UMAP ##

p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = ovgenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)
######### Marker Genes Imputation with MAGIC
brain2 <- addImputeWeights(brain2)
p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = ovgenes, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(brain2)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain2, 
  groupBy = "Clusters", 
  geneSymbol = ovgenes, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

##################################################################
######### markerGenes.mic UMAP ##
##################################################################
markerGenes.mic
markerGenes.mic=markerGenes.mic[-2]
p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)
######### Marker Genes Imputation with MAGIC
p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(brain2)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes.mic, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Mic_Marker-Genes.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

#saveArchRProject(ArchRProj = brain, outputDirectory = "Save-brain", load = FALSE)
saveRDS(brain2,file='brain.microglia.filter/brain2.rds')

##################################################################
######### cross-platform linkage of scATAC-seq cells with scRNA-seq cells
################################################################## 

brain.rna=readRDS('ROSMAP.Microglia.6regions.seurat.harmony.selected.rds')
table(brain.rna@meta.data$seurat_clusters)
brain.rna$groupRNA=paste('MG',brain.rna@meta.data$seurat_clusters,sep='')

brain3 <- addGeneIntegrationMatrix(
  ArchRProj = brain2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = brain.rna,
  addToArrow = F,
  groupRNA = "groupRNA",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(brain3$Clusters, brain3$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))
table(brain3$Clusters)

unique(brain3$predictedGroup_Un)

p1 <- plotEmbedding(
  brain3, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
)

plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = brain3, addDOC = FALSE, width = 5, height = 5)

getAvailableMatrices(brain3)

##### 
brain3 <- addImputeWeights(brain3)
meta=brain3@cellColData
write.table(meta,file='brain.microglia.filter/brain3.meta.txt',sep='\t',quote=F)
saveRDS(brain3,file='brain.microglia.filter/brain3.rds')

table(meta$Clusters)
remove<-c("C2","C3") ## remove clusters with very few cells
keep<-!(brain3$Clusters %in% remove)
sum(keep=="TRUE") ## the number of cell to be kept
brain4=brain3[keep,] ## remove small clusters

remapClust <- c("C1" = "C1","C4" = "C2","C5" = "C3","C6" = "C4","C7" = "C5","C8" = "C6","C9" = "C7","C10" = "C8","C11" = "C9")
labelNew <- mapLabels(names(remapClust), oldLabels = names(remapClust), newLabels = remapClust)
labelNew

brain4$Clusters <- mapLabels(brain4$Clusters, newLabels = labelNew, oldLabels = names(remapClust))


p1 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "region", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "projid", embedding = "UMAP")

plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters.filtered.pdf", ArchRProj = brain4, addDOC = FALSE, width = 6, height = 6)


#### plot marker genes 

##### identify marker genes ############################################
##############################################################
### 
markersGS <- getMarkerFeatures(
  ArchRProj = brain4, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList,file='brain.microglia.filter/markerList.rds')

df=read.table('microglia.markers.human.txt',sep = '\t')
markerGenes.mic = as.character(df$V1)
rnadeg=read.table('ROSMAP.Microglia.6regions.seurat.harmony.clusterDEGs.txt',header = T,sep = '\t')
markerGenes=unique(as.character(rnadeg$gene))

######
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)

##
heatmapGS.mic <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
  labelMarkers = markerGenes.mic,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)
############
heatmapGS.both <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = c(markerGenes,markerGenes.mic),
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS.both, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker-Heatmap", width = 20, height = 6, ArchRProj = brain2, addDOC = FALSE)

## overlap between snRNA subtype markers and snATAC subtype markers
atacgenes=colnames(heatmapGS@matrix)
ovgenes=intersect(atacgenes,markerGenes)

##################################################################
######### markerGenes.mic UMAP ##
##################################################################
markerGenes.mic
markerGenes.mic=markerGenes.mic[-2]
mymat=getMatrixFromProject(
  ArchRProj = brain4,
  useMatrix = "GeneScoreMatrix")
brain4 <- addImputeWeights(brain4)
matGS <- imputeMatrix(assay(mymat), getImputeWeights(brain4))
dgc_imput_mat <- as(matGS, "dgCMatrix")
str(dgc_imput_mat)
rownames(dgc_imput_mat)=mymat@elementMetadata$name
saveRDS(dgc_imput_mat,file='brain.microglia.filter/brain4.imput_genescoremat.rds')

umaps=getEmbedding(ArchRProj = brain4, embedding = "UMAP", returnDF = TRUE)
saveRDS(umaps,file='brain.microglia.filter/brain4.umaps.rds')


###################

p <- plotEmbedding(
  ArchRProj = brain4, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)
######### Marker Genes Imputation with MAGIC
p <- plotEmbedding(
  ArchRProj = brain4, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(brain4)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain4, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes.mic, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Mic_Marker-Genes.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)

##########################################################
## 

meta=brain4@cellColData
write.table(meta,file='brain.microglia.filter/brain4.meta.txt',sep='\t',quote=F)
saveRDS(brain4,'brain.microglia.filter/brain4.rds')


#### 08/15/2022 #### add cell state
brain=readRDS('brain.microglia.filter/brain4.rds')
getAvailableMatrices(brain)
ids=c('homeostatic','activated1','activated1','homeostatic','homeostatic','activated2','homeostatic','activated2','activated2')
names(ids)=paste('C',1:9,sep='')
brain$cellstate=ids[brain$Clusters]

brain <- addPeak2GeneLinks(ArchRProj = brain,
  useMatrix = "GeneScoreMatrix", ## by default is GeneIntegrationMatrix
  reducedDims = "IterativeLSI")


### pcc as 0.3 to get links
p2g <- getPeak2GeneLinks(
    ArchRProj = brain,
    corCutOff = 0.3,
    resolution = 1,
    returnLoops = F
)
p2g[[1]]

p <- plotPeak2GeneHeatmap(ArchRProj = brain, corCutOff = 0.3,k=6,groupBy = "cellstate",returnMatrices = T)
plotPDF(p, 
    name = "plotPeak2GeneHeatmap_cellstate.pdf", 
    ArchRProj = brain, 
    addDOC = FALSE, width = 12, height = 12)
saveRDS(p,file='brain.microglia.filter/Peak2GeneLinks.heatmap.matrix.rds')


df=read.table('microglia.markers.human.txt',sep = '\t')
markerGenes.mic = as.character(df$V1)  

p <- plotBrowserTrack(
    ArchRProj = brain, 
    groupBy = "cellstate", 
    geneSymbol = markerGenes.mic, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(
    ArchRProj = brain,
    corCutOff = 0.3,
    resolution = 1)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = brain, 
    addDOC = FALSE, width = 5, height = 5)
### call peaks based on new annotation
pathToMacs2 <- findMacs2()
## Making Pseudo-bulk Replicates
brain <- addGroupCoverages(ArchRProj = brain, groupBy = "cellstate")
brain <- addReproduciblePeakSet(ArchRProj = brain, groupBy = "cellstate", pathToMacs2 = pathToMacs2)
getPeakSet(brain)

saveRDS(getPeakSet(brain),file='brain.microglia.filter/getPeakSet_cellstate.rds')

## Add Peak Matrix

brain <- addPeakMatrix(brain)
getAvailableMatrices(brain)
## marker peaks
markersPeaks <- getMarkerFeatures(
    ArchRProj = brain, 
    useMatrix = "PeakMatrix", 
    groupBy = "cellstate",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
saveRDS(brain,'brain.microglia.filter/microglia.snATAC.ArchRobj.rds')
markersPeaks
saveRDS(markersPeaks,file='brain.microglia.filter/markersPeaks_cellstate.rds')

markerList <- getMarkers(markersPeaks, cutOff = "Pval <= 0.05 & Log2FC >= 0.5")
markerList
str(markerList)
saveRDS(markerList,file='brain.microglia.filter/markersPeaks_cellstate.pval05.rds')

markerList <- getMarkers(markersPeaks, cutOff = "Pval <= 0.01 & Log2FC >= 0.5")
markerList
str(markerList)
saveRDS(markerList,file='brain.microglia.filter/markersPeaks_cellstate.pval01.rds')

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "Pval <= 0.01 & Log2FC >= 0.5",
  transpose = FALSE
)

plotPDF(heatmapPeaks, name = "Peak-cellstate_Marker-Heatmap.pval01", width = 8, height = 6, ArchRProj = brain, addDOC = FALSE)


heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5",
  transpose = FALSE
)

plotPDF(heatmapPeaks, name = "Peak-cellstate_Marker-Heatmap.pval05", width = 8, height = 6, ArchRProj = brain, addDOC = FALSE)

##### gene markers based on gene score

##### identify marker genes ############################################
##############################################################
### 
markersGS <- getMarkerFeatures(
  ArchRProj = brain, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "cellstate",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "Pval <= 0.01 & Log2FC >= 0.25")

saveRDS(markerList,file='brain.microglia.filter/markerList_cellsate.rds')

df=read.table('microglia.markers.human.txt',sep = '\t')
markerGenes.mic = as.character(df$V1)
rnadeg=read.table('ROSMAP.Microglia.6regions.seurat.harmony.clusterDEGs.txt',header = T,sep = '\t')
markerGenes=unique(as.character(rnadeg$gene))

######
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "Pval <= 0.01 & Log2FC >= 0.25", 
  labelMarkers = markerGenes,
  transpose = FALSE
)

#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker_cellstate-Heatmap", width = 12, height = 6, ArchRProj = brain, addDOC = FALSE)

##
heatmapGS.mic <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "Pval <= 0.01 & Log2FC >= 0.25", 
  labelMarkers = markerGenes.mic,
  transpose = FALSE
)

#ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker_cellstate-Heatmap", width = 12, height = 6, ArchRProj = brain, addDOC = FALSE)
############
heatmapGS.both <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "Pval <= 0.01 & Log2FC >= 0.25", 
  labelMarkers = c(markerGenes,markerGenes.mic),
  transpose = FALSE
)

ComplexHeatmap::draw(heatmapGS.both, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker_cellstate-Heatmap", width = 20, height = 6, ArchRProj = brain, addDOC = FALSE)


######### 

