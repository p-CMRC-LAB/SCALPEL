

#Libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(argparse)
library(scales)



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
dtab = reads %>% distinct(frag_id,read_id,dist_END) %>% arrange(dist_END)
dtab = table(dtab$dist_END) %>% data.table()
dtab$V1 = as.numeric(dtab$V1)
colnames(dtab) = c('transcriptomic_distance','counts')

#let's plot the fragments ends positions and reads
# ggplot(dtab, aes(transcriptomic_distance, counts)) +
#     geom_line() + geom_area(fill="cornflowerblue") + geom_smooth(span=0.05, color="red")

#let's set interval axis
part_neg = c(seq(0,dtab$transcriptomic_distance[1],-BINS),dtab$transcriptomic_distance[1]) %>% unique() %>% rev()
part_pos = c(seq(0, max(dtab$transcriptomic_distance), BINS), max(dtab$transcriptomic_distance)) %>% unique()
intervals = unique(c(part_neg,part_pos))
intervals_counts = list()
intervals_idx = list()
intervals_probs = list()
for(i in 1:(length(intervals)-1)){
    left_born = intervals[i]
    right_born = intervals[i+1]
    #index
    intervals_idx[[i]] = seq(left_born, right_born)
    #counts
    intervals_counts[[i]] = rep((dtab %>% filter(transcriptomic_distance>= left_born & transcriptomic_distance < right_born))$counts %>% sum(), length(intervals_idx[[i]]))
}
interval.tab = data.table(transcriptomic_distance = unlist(intervals_idx), counts_cum = unlist(intervals_counts)) %>% distinct(transcriptomic_distance,.keep_all = T)
interval.tab$probs_bin = (interval.tab$counts_cum / sum(interval.tab$counts_cum)) * 100
interval.tab = left_join(interval.tab, dtab)
interval.tab$counts = interval.tab$counts %>% tidyr::replace_na(0)
#let's normalize for visualization
interval.tab$probability_scaled = interval.tab$probs_bin %>% scales::rescale(to = c(0,max(interval.tab$counts)))
#interval.tab$loess = loess(counts~transcriptomic_distance, data = interval.tab, span = 0.1, weights = interval.tab$probs_bin)$fitted %>% scales::rescale(to=c(0,max(interval.tab$counts)))

ggplot(interval.tab) +
    geom_area(aes(transcriptomic_distance,counts), fill="cornflowerblue", size=0.5) +
#    geom_line(aes(transcriptomic_distance,loess), color="gray20", color="red", size=1) +
    geom_line(aes(transcriptomic_distance,probability_scaled), color="red",size=1) +
    theme_classic() +
    ggtitle("Fragments distribution on transcriptomic space")
ggsave(args$OUTPUT_PDF, scale = 1, device = "pdf",units = "in", width = 11.69, height = 8.27)


# [Writing]
# --------
#interval.tab$probs_bin = interval.tab$loess
ggsave(args$OUTPUT_PDF, scale = 1, device = "pdf",units = "in", width = 11.69, height = 8.27)
fwrite(interval.tab[,c('transcriptomic_distance','counts','probs_bin')], file = args$OUTPUT_PROB, sep="\t", col.names = F)
