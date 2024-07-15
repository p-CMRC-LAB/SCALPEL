library(dplyr)
library(data.table)
library(argparse)
library(stringr)
library(tidyr)


#Argument Parser
#===============
parser = ArgumentParser(description='Merging DGE & isoforms Prediction')
parser$add_argument('DGE', type="character", help='path of DGE file')
parser$add_argument('PRED', type="character", help='path of predicted isoforms')
parser$add_argument('OUTP', type="character", help="path of output file")
args = parser$parse_args()

#Files opening
preds = fread(args$PRED, col.names=c("bc","gene_name","transcript_name","tr_prob"))
dge = fread(args$DGE)

#melting
dge = melt(dge)

#renaming
colnames(dge) = c("gene_name","bc","counts")

#left_joining by barcodes & gene_name
dge_preds = left_join(preds, dge) %>% stats::na.omit() %>% mutate(tr_counts=tr_prob * counts, gene_transcript = paste0(gene_name, "***", transcript_name)) %>% distinct(gene_transcript, bc, tr_counts) %>% pivot_wider(names_from=bc, values_from=tr_counts, values_fill = 0)


#writing
fwrite(dge_preds, file=args$OUTP, sep="\t")
