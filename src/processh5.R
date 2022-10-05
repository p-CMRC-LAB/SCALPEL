

#libraries
#*********
library(Seurat)
library(argparse)
library(data.table)

#settings
#********
options(future.globals.maxSize = 1500 * 1024^2)


#Argunent Parser
#**************
parser = ArgumentParser(description='Process h5 file 10X... ')
parser$add_argument('h5file', metavar='PRFfip', type="character",
                    help='path of h5 file')
parser$add_argument('output', metavar='PPFfop', type="character",
                    help='path of output file')
args = parser$parse_args()

print("file opening...")
df = as.data.frame(as.matrix(Read10X_h5(args$h5file)))
df = cbind(data.frame(GENE = rownames(df)), df)


#writing
fwrite(df, file = args$output, sep="\t", row.names=FALSE)
