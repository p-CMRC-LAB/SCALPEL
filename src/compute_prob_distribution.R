

#Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(argparse)
library(scales)
# library(distributionsrd)
# library(fitdistrplus)



#Argunent Parser
#**************
parser = ArgumentParser(description='Probability of unique transcript')
parser$add_argument('PATH_OF_UNIQUE_TR', metavar='Tfip', type="character", help='path of unique isoform all reads')
parser$add_argument('GFRAC', metavar='gfrac', type="character", help='fraction of gene to keep in filtering')
parser$add_argument('BINS', metavar='bins', type="character", help='size of bins')
parser$add_argument('OUTPUT_PROB', metavar='Tfop', type="character", help='path of output unique isoform intervals probabilities')
parser$add_argument('OUTPUT_PDF', metavar='PDFfop', type="character", help='path of output pdf pictures')
args = parser$parse_args()

#Variables
QUANTILE_THRESHOLD = args$GFRAC
BINS = as.numeric(args$BINS)

#Files Opening (1)
#*************
#get list of of unique read file
files = list.files(args$PATH_OF_UNIQUE_TR, pattern = "*.fragment_filtered_unique", full.names = T)
reads = parallel::mclapply(files, function(x) fread(x), mc.preschedule = T, mc.cores = 5)
reads = do.call(rbind, reads)
reads

#filter out genes bringing variabilities
genes_counts = table(reads$gene_name)
genes_quants = genes_counts %>% quantile(seq(0,1,0.01))
reads = reads %>% filter(gene_name %in% names(genes_counts[genes_counts<genes_quants[[QUANTILE_THRESHOLD]]]))

#get distinct readid / 3end
dtab = reads %>% distinct(read_id,dist_END) %>% arrange(dist_END)
dtab = table(dtab$dist_END) %>% data.table()
dtab$V1 = as.numeric(dtab$V1)
colnames(dtab) = c('transcriptomic_distance','counts')

#create a table covering the whole transcriptomic space
MIN = min(reads$dist_END)
MAX = max(reads$dist_END)
dtab_cov = data.table(data.frame(transcriptomic_distance = seq(MIN,MAX)))
#merge
merged = left_join(dtab_cov, dtab)
merged$counts = merged$counts %>% tidyr::replace_na(0)
#eliminate duplicates
merged = merged %>% distinct(.keep_all=TRUE)

#plot
# ggplot(merged, aes(transcriptomic_distance,counts)) +
#   geom_point(size=1) +
#   geom_line(size=0.5) +
#   geom_smooth(span=0.30) +
#   theme_bw()

#breaks
if(MIN < 0){
	breaks = unique(c(-Inf,sort(unique(c(seq(0,MIN,-BINS),seq(0,MAX,BINS)))),MAX)); breaks
}else{
	breaks = unique(c(-Inf,sort(unique(c(seq(0,MAX,BINS)))),MAX)); breaks
}

#get probability distribution on BINS intervals
merged.2 = do.call(rbind, lapply(1:(length(breaks)-1), function(x){
  A = breaks[x]
  B = breaks[x+1]
  ends = merged$transcriptomic_distance[(merged$transcriptomic_distance>A) & (merged$transcriptomic_distance<=B)]
  if(length(ends)==0){ends = B}
  ends_counts = merged$counts[merged$transcriptomic_distance %in% ends]
  if(length(ends_counts)==0){ends = 0}
  return(data.frame(transcriptomic_distance = ends,bin_counts = sum(ends_counts)))
}))
merged.2 = merged.2 %>% distinct(transcriptomic_distance, .keep_all=TRUE)

#final_merging
merged.final = left_join(merged.2, merged) %>% data.table(); merged.final

#final scaling
merged.final$bin_counts_pb = merged.final$bin_counts / sum(merged.final$bin_counts)
merged.final$probability_normalized = merged.final$bin_counts_pb %>% scales::rescale(to=c(0,1))
merged.final$counts_normalized = merged.final$counts %>% scales::rescale(to=c(0,1))
merged.final = merged.final[,c('transcriptomic_distance','counts','bin_counts_pb','counts_normalized','probability_normalized','bin_counts')]
colnames(merged.final) = c('transcriptomic_distance','counts','probs_bin','counts_normalized','probability_normalized','bin_counts')
merged.final$probs_bin = merged.final$probs_bin * 100

merged.loess = loess(probs_bin~transcriptomic_distance, data = merged.final[merged.final$transcriptomic_distance>=0], span = 0.05)
merged.final$loess = 0
merged.final$loess[merged.final$transcriptomic_distance>=0] = merged.loess$fitted


#plot
ggplot(merged.final) +
  geom_line(aes(transcriptomic_distance,probability_normalized, col='bins distribution')) +
  geom_line(aes(transcriptomic_distance,loess, col='loess fitting of bins distribution')) +
  geom_line(aes(transcriptomic_distance,counts_normalized, col='read counts')) +
  theme_bw(base_size = 15)


#writing
merged.final$loess = merged.final$loess * 100
merged.final$probs_bin = merged.final$loess
ggsave(args$OUTPUT_PDF, scale = 2, device = "jpeg")
fwrite(merged.final[,c('transcriptomic_distance','counts','probs_bin','counts_normalized','probability_normalized')], file = args$OUTPUT_PROB, sep="\t")
