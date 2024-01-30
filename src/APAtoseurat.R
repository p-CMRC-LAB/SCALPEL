

library(argparse)
library(data.table)
library(dplyr)
library(Seurat)


#ArgumentParser
#==============
parser = ArgumentParser(description='Process counts file')
parser$add_argument('APA_PATH', type="character", help='path of APA counts file')
parser$add_argument('sample_name', type="character", help='sample_name')
parser$add_argument('OUTPUT', type="character", help='path of output seurat file')
args = parser$parse_args()

#file opening and writing
read.table(args$APA_PATH, sep="\t", header = T, row.names = 1) %>% 
    Seurat::CreateSeuratObject(names.field = 1, min.cells = 10, min.features = 4, project = args$sample_name) %>%
    saveRDS(file=paste0(args$sample_name,"_seurat.RDS"))

