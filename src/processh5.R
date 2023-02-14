

#libraries
#*********
# library(Seurat)
library(argparse)
library(data.table)
library(Matrix)

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

readh5 = function (filename, use.names = TRUE, unique.features = TRUE)
{
    if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Please install hdf5r to read HDF5 files")
    }
    if (!file.exists(filename)) {
        stop("File not found")
    }
    infile <- hdf5r::H5File$new(filename = filename, mode = "r")
    genomes <- names(x = infile)
    output <- list()
    if (hdf5r::existsGroup(infile, "matrix")) {
        if (use.names) {
            feature_slot <- "features/name"
        }
        else {
            feature_slot <- "features/id"
        }
    }
    else {
        if (use.names) {
            feature_slot <- "gene_names"
        }
        else {
            feature_slot <- "genes"
        }
    }
    for (genome in genomes) {
        counts <- infile[[paste0(genome, "/data")]]
        indices <- infile[[paste0(genome, "/indices")]]
        indptr <- infile[[paste0(genome, "/indptr")]]
        shp <- infile[[paste0(genome, "/shape")]]
        features <- infile[[paste0(genome, "/", feature_slot)]][]
        barcodes <- infile[[paste0(genome, "/barcodes")]]
        sparse.mat <- sparseMatrix(i = indices[] + 1, p = indptr[],
            x = as.numeric(x = counts[]), dims = shp[], giveCsparse = FALSE)
        if (unique.features) {
            features <- make.unique(names = features)
        }
        rownames(x = sparse.mat) <- features
        colnames(x = sparse.mat) <- barcodes[]
        sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
        if (infile$exists(name = paste0(genome, "/features"))) {
            types <- infile[[paste0(genome, "/features/feature_type")]][]
            types.unique <- unique(x = types)
            if (length(x = types.unique) > 1) {
                message("Genome ", genome, " has multiple modalities, returning a list of matrices for this genome")
                sparse.mat <- sapply(X = types.unique, FUN = function(x) {
                  return(sparse.mat[which(x = types == x), ])
                }, simplify = FALSE, USE.NAMES = TRUE)
            }
        }
        output[[genome]] <- sparse.mat
    }
    infile$close_all()
    if (length(x = output) == 1) {
        return(output[[genome]])
    }
    else {
        return(output)
    }
}

print("file opening...")
df = as.data.frame(as.matrix(readh5(args$h5file)))
df = cbind(data.frame(GENE = rownames(df)), df)


#writing
fwrite(df, file = args$output, sep="\t", row.names=FALSE)
