

library(argparse)
library(data.table)
library(dplyr)
library(Seurat)


#ArgumentParser
#==============
parser = ArgumentParser(description='Process 10X repo file')
parser$add_argument('repo_path', type="character", help='path of 10x repo')
args = parser$parse_args()


#define Sample_name
sample_name = basename(args$repo_path)

#Reading
s.mat = data.frame(Read10X(paste0(args$repo_path,"/outs/filtered_feature_bc_matrix/")))


#Writing
#Write Barcodes file
barcodes.tab = data.table(colnames(s.mat)) %>%
  mutate(V1=stringr::str_replace(V1,"\\.","\\-"))
barcodes.tab %>%
  fwrite(file = paste0(sample_name, ".barcodes.txt"), col.names = F)

#write Matrix
colnames(s.mat) = barcodes.tab$V1 
s.mat = s.mat %>% data.table(keep.rownames=T)
colnames(s.mat)[1] = "GENE"
s.mat %>% fwrite(file = paste0(sample_name, ".counts.txt"), sep="\t", row.names=F)