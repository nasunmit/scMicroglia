
###
load('ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.degset.RData')
str(degset)
### 
library(enrichR)
options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
dbs <- listEnrichrDbs()

cutoff=0.05
mydbs <- c("TRANSFAC_and_JASPAR_PWMs","ChEA_2022",'ENCODE_TF_ChIP-seq_2015')
for (c in names(degset)){
  genes=as.character(degset[[c]])
  hspgenes=c(genes[grep('^HSP',genes)],genes[grep('^RPL',genes)],genes[grep('^RPS',genes)])
  genes=setdiff(genes,hspgenes)
  if (length(genes)>89){genes=genes[1:89]}
  enriched <- enrichr(genes, mydbs)
  
  out=enriched$TRANSFAC_and_JASPAR_PWMs
  sel=out[out$Adjusted.P.value<cutoff,]

  out=enriched$ChEA_2022
  sel=rbind(sel,out[out$Adjusted.P.value<cutoff,])
  
  out=enriched[['ENCODE_TF_ChIP-seq_2015']]
  sel=rbind(sel,out[out$Adjusted.P.value<cutoff,])
  
  write.table(sel,file=paste(paste('microglia.subtype.markers',c,sep='.'),'.enrichedTFs.txt',sep=''),sep='\t',quote = F)
  
}

