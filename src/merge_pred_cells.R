

#script for generating from the prob files, the corresponding cell files


#import libs
library(data.table)
library(dplyr)
library(argparse)


#Argunent Parser
#**************
parser = ArgumentParser(description='Merge predictions...')
parser$add_argument('PATH_OF_READS_FILE', metavar='PRFfip', type="character",
                    help='path of predicted fragments')
parser$add_argument('PATH_OF_MERGED_FILE', metavar='PMFfop', type="character",
                    help='path of predicted fragments merged')
parser$add_argument('nTHREADS', metavar='nThreads', type="character",
                    help='path of predicted fragments merged')
args = parser$parse_args()

#files import
pred_files <- list.files(args$PATH_OF_READS_FILE, full.names = T, pattern = "*.pred")

#Files Opening
print("Files reading...")
RES_all <- lapply(pred_files, function(x){
  print(x)
  #reading
  pred = fread(x, nThread =1)
}) %>% rbindlist()


#writing
print("Writing...")
fwrite(RES_all, file = args$PATH_OF_MERGED_FILE, quote = F, col.names = T, row.names = F, sep = "\t", nThread = as.numeric(args$nTHREADS))
