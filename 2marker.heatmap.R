##########################
library(scales)
library(plyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
wr <- colorRampPalette(colors = c( "white", "red"))(100)
rwb <- colorRampPalette(colors = c("blue", "white", "red"))(100)
ryb=colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)
prgn=colorRampPalette(rev(brewer.pal(n = 10, name ="PRGn")))(100)
wb <- colorRampPalette(colors = c( "white", "blue"))(100)
##
id='ROSMAP.Microglia.6regions.seurat.harmony.selected'
######## DEGs 
load('ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.degset.RData')
degs=read.table(file='ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt',sep = '\t',header=T)
clusters=paste('MG',degs$cluster,sep='')
clusters=factor(clusters,levels = unique(clusters))

bulkdata=read.table('ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs_pseudobulk.txt',sep = '\t',header = T)
colnames(bulkdata)=gsub('C','MG',colnames(bulkdata))
head(bulkdata)

data=c()
nums=c()
for (i in names(degset)){
  data=rbind(data,bulkdata[degset[[i]],])
  nums=c(nums,length(degset[[i]]))
}

library(ComplexHeatmap)
library(circlize)
mycol=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB")
test=t(scale(t(data)))
leftanno = rowAnnotation(size = anno_block(gp = gpar(fill = mycol),labels_rot=0,
                                           labels = as.character(nums), 
                                           labels_gp = gpar(col = "white", fontsize = 10)))
ht=Heatmap(scale(test),cluster_columns = F,name = 'z-score',row_split = clusters,cluster_rows = F,
           row_title_rot = 0,left_annotation=leftanno,border = TRUE,border_gp = gpar(col = "black", lwd = 0.5),
           show_row_names = F,row_title_gp = gpar(fontsize = 8))

pdf(file='ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.complexheatmap.withNum.pdf',height = 10,width = 8)
ht
dev.off()





