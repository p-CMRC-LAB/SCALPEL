

library(argparse)
library(dplyr)
library(doParallel)
library(data.table)

#ArgumentParser
#==============
parser = ArgumentParser(description='Process GTF file')
parser$add_argument('GTF_PATH', type="character", help='path of gtf file')
args = parser$parse_args()

#File opening
#============
print("file opening...")
gtf = rtracklayer::import(args$GTF_PATH) %>% data.frame() %>% filter(type=="exon")

#writing
#=======
print("file writing...")
res = lapply(split(gtf, gtf$seqnames), function(x){
    fwrite(x, file = paste0(x$seqnames[1],".gtf"), sep="\t")
})
