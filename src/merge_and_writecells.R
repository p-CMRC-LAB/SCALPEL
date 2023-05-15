



library(data.table)
library(dplyr)
library(doParallel)
library(argparse)



#Argunent Parser
#**************
parser = ArgumentParser(description='merge cells from chromosome files')
parser$add_argument('frags_path', type="character", help='path of unique isoform all reads')
args = parser$parse_args()


#get files
all_files = list.files(args$frags_path, pattern = "*.frag_prob", full.names = T)

#run awk on each of those files
command.awk = lapply(all_files, function(x) paste0("awk '{print >> $1\".cell\"}' ", x)) %>% unlist()
res = lapply(command.awk, function(x){
  system(x)
})






