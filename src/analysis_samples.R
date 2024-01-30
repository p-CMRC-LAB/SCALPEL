

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


#get RDS files
all_files = list.files(path=args$SEURAT_OBJ_LOC, pattern = "*.RDS", full.names = T)
print(all_files)

#Opening (1)
#*******
all.objs = lapply(all_files, function(x){
  print(x)
  curr.obj = readRDS(x)
})

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
	cluster_df = fread(args$CLUSTERS, col.names = c("Barcodes", "Clusters"))

	#b merge cluster into the seurat object
	seurat.obj$Clusters = (dplyr::left_join(seurat.obj@meta.data, cluster_df) %>% distinct(Barcodes, .keep_all=T))$Clusters
}

#Differential analysis (3)
#*********************
my_isoforms = Find_isoforms(seurat.obj, condition="Clusters", pval_adjusted=0.05, threshold_tr_abundance = 0.20)


#Writing differential analysis table (4)
#***********************************
fwrite(my_isoforms, file=args$OUTPUT, sep="\t")
saveRDS(seurat.obj, file="final_seurat_obj.RDS")