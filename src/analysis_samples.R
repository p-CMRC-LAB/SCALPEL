

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


<<<<<<< HEAD
THR_ABUND = 0.10
THR_PVAL = 0.05


#get RDS files
all_files = list.files(path=args$SEURAT_OBJ_LOC, pattern = "*.RDS", full.names = T)
=======
#get RDS files
all_files = list.files(path=args$SEURAT_OBJ_LOC, pattern = "*.RDS", full.names = T)
print(all_files)
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

#Opening (1)
#*******
all.objs = lapply(all_files, function(x){
  print(x)
  curr.obj = readRDS(x)
})

<<<<<<< HEAD
print(all.objs)
=======
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
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
<<<<<<< HEAD
	cluster_df = fread(args$CLUSTERS, col.names = c("orig.ident", "Barcodes", "Clusters"))
=======
	cluster_df = fread(args$CLUSTERS, col.names = c("Barcodes", "Clusters"))
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

	#b merge cluster into the seurat object
	seurat.obj$Clusters = (dplyr::left_join(seurat.obj@meta.data, cluster_df) %>% distinct(Barcodes, .keep_all=T))$Clusters
}

<<<<<<< HEAD
print(seurat.obj)
print( seurat.obj[[]] %>% head() )
print( seurat.obj[[]] %>% tail() )

=======
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
#Differential analysis (3)
#*********************
if(n_distinct(seurat.obj$Clusters) == 1){
	my_isoforms = data.table(output="Provide more condition with levels > 1 for differential isoform usage analysis")
}else{
<<<<<<< HEAD
        print("in")
	my_isoforms = Find_isoforms(seurat.obj, condition="Clusters", pval_adjusted=THR_PVAL, threshold_tr_abundance=THR_ABUND)
=======
	my_isoforms = Find_isoforms(seurat.obj, condition="Clusters", pval_adjusted=0.05, threshold_tr_abundance = 0.20)
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
}

#Writing differential analysis table (4)
#***********************************
fwrite(my_isoforms, file=args$OUTPUT, sep="\t", nThread=1)
<<<<<<< HEAD
saveRDS(seurat.obj, file="iDGE_seurat.RDS")
=======
saveRDS(seurat.obj, file="final_seurat_obj.RDS")
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
