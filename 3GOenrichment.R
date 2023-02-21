#####

##########
df=read.table('ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt',header = T,sep = '\t')
clusters=names(table(df$cluster))
degset=list()
for (c in clusters){
  genes=as.character(df[df$cluster==c,7])
  degset[[paste('MG',c,sep='')]]=genes
}
str(degset)

#### go enrichment 
### 
library(enrichR)
options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
dbs <- listEnrichrDbs()

dbs[grep('KEGG',dbs$libraryName),]
dbs[grep('GO_Biological_Process',dbs$libraryName),]

cutoff=0.05
mydbs <- c("GO_Biological_Process_2021","KEGG_2021_Human")
for (c in names(degset)){
  genes=as.character(degset[[c]])
  enriched <- enrichr(genes, mydbs)
  
  out=enriched$GO_Biological_Process_2021
  out=out[out$Adjusted.P.value<cutoff,]
  write.table(out,file=paste(c,'.enrichedGO.txt',sep=''),sep='\t',quote = F)
  
  out=enriched$KEGG_2021_Human
  out=out[out$Adjusted.P.value<cutoff,]
  write.table(out,file=paste(c,'.enrichedKEGG.txt',sep=''),sep='\t',quote = F)
  
}

