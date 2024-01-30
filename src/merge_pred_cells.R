

#script for generating from the prob files, the corresponding cell files


#import libs
library(data.table)
library(dplyr)
library(parallel)
library(argparse)


#Argunent Parser
#**************
parser = ArgumentParser(description='Merge predictions...')
parser$add_argument('PATH_OF_READS_FILE', metavar='PRFfip', type="character",
                    help='path of predicted fragments')
parser$add_argument('PATH_OF_MERGED_FILE', metavar='PMFfop', type="character",
                    help='path of predicted fragments merged')
args = parser$parse_args()

#files import
pred_files <- list.files(args$PATH_OF_READS_FILE, full.names = T, pattern = "*.pred")

#Files Opening
print("Files reading...")
RES <- mclapply(pred_files, function(x){
  print(x)
  #reading
  pred = fread(x)
}, mc.preschedule = T, mc.cores = 1)


#Merging
print("Files merging")
RES_all = do.call(rbind, RES)


#writing
print("Writing...")
fwrite(RES_all, file = args$PATH_OF_MERGED_FILE, quote = F, col.names = T, row.names = F, sep = "\t", nThread = 20)
