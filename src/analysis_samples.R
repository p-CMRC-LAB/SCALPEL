

library(argparse)
library(data.table)
suppressWarnings(library(dplyr))
library(Seurat) %>% suppressWarnings()
source("scalpel_library.R") %>% suppressPackageStartupMessages()



#ArgumentParser
#==============
parser = ArgumentParser(description='Process seurat file')
parser$add_argument('SEURAT_OBJ_LOC', type="character", help='path of APA counts file')
parser$add_argument('CLUSTERS', type="character", help='path of clusters file')
parser$add_argument('OUTPUT', type="character", help='path of output seurat file')
args = parser$parse_args()


THR_ABUND = 0.10
THR_PVAL = 0.05


#get RDS files
all_files = list.files(path=args$SEURAT_OBJ_LOC, pattern = "*.RDS", full.names = T)

#Opening (1)
#*******
all.objs = lapply(all_files, function(x){
  print(x)
  curr.obj = readRDS(x)
})

print(all.objs)
#merging
if(length(all.objs)>1){
	seurat.obj = merge(all.objs[[1]], all.objs[2:length(all.objs)])
}else{
	seurat.obj = all.objs[[1]]
}

#annotate cell barcodes with clusters if needed (default samples) (2)
#****************************************************************
seurat.obj$Barcodes = rownames(seurat.obj@meta.data)

if(args$CLUSTERS == "NULL"){
	seurat.obj$Clusters = seurat.obj$orig.ident
}else{
	#a. read clusters table
	cluster_df = fread(args$CLUSTERS, col.names = c("orig.ident", "Barcodes", "Clusters"))

	#b merge cluster into the seurat object
	seurat.obj$Clusters = (dplyr::left_join(seurat.obj@meta.data, cluster_df) %>% distinct(Barcodes, .keep_all=T))$Clusters
}

print(seurat.obj)
print( seurat.obj[[]] %>% head() )
print( seurat.obj[[]] %>% tail() )

#Differential analysis (3)
#*********************
if(n_distinct(seurat.obj$Clusters) == 1){
	my_isoforms = data.table(output="Provide more condition with levels > 1 for differential isoform usage analysis")
}else{
        print("in")
	my_isoforms = Find_isoforms(seurat.obj, condition="Clusters", pval_adjusted=THR_PVAL, threshold_tr_abundance=THR_ABUND)
}

#Writing differential analysis table (4)
#***********************************
fwrite(my_isoforms, file=args$OUTPUT, sep="\t", nThread=1)
saveRDS(seurat.obj, file="iDGE_seurat.RDS")
