
suppressMessages(suppressWarnings(library(scales)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(ggplot2)))



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


#1. Files Opening
#++++++++++++++++
print("opening...")
#get list of of unique read file
reads = fread(args$PATH_OF_UNIQUE_TR, col.names=c("seqnames.rd","start.rd","end.rd","strand.rd","dist_END","frag.id","start","end","gene_name","transcript_name","bulk_TPMperc"), nThread = 1)

#2. Processing
#+++++++++++++
print("processing...")
#- get gene counts
gene_counts = reads %>%
  distinct(seqnames.rd,start.rd,end.rd,strand.rd,frag.id,gene_name,transcript_name) %>%
  group_by(gene_name) %>%
  summarise(nb.reads = n()) %>%
  data.table()
gene.quantiles = quantile(gene_counts$nb.reads, seq(0,1,0.01))
gene.tokeep = dplyr::filter(gene_counts, nb.reads < gene.quantiles[[QUANTILE_THRESHOLD]])$gene_name
#filtering
reads = reads %>%
  dplyr::filter(gene_name %in% gene.tokeep) %>%
  mutate(dist_END = as.numeric(dist_END))

#- get the number of distinct reads at each 3'end position
read_tab = distinct(reads, dist_END, start.rd, end.rd) %>%
  group_by(dist_END) %>%
  reframe(read.counts = n()) %>%
  arrange(dist_END) %>%
  na.omit() %>%
  data.table()
colnames(read_tab) = c("transcriptomic_distance", "counts")

#let's plot the reads ends positions and reads
ggplot(read_tab, aes(transcriptomic_distance, counts)) +
  geom_line() + geom_area(fill="cornflowerblue") + geom_smooth(span=0.05, color="red")

#let's set interval axis
part_neg = c(seq(0,read_tab$transcriptomic_distance[1],-BINS),read_tab$transcriptomic_distance[1]) %>% unique() %>% rev()
part_pos = c(seq(0, max(read_tab$transcriptomic_distance), BINS), max(read_tab$transcriptomic_distance)) %>% unique()
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
  intervals_counts[[i]] = rep(
    (read_tab %>% dplyr::filter(transcriptomic_distance>= left_born & transcriptomic_distance < right_born))$counts %>%
      sum(), length(intervals_idx[[i]]))
}
interval.tab = data.table(transcriptomic_distance = unlist(intervals_idx), counts_cum = unlist(intervals_counts)) %>%
  distinct(transcriptomic_distance,.keep_all = T)
interval.tab$probs_bin = (interval.tab$counts_cum / sum(interval.tab$counts_cum)) * 100
interval.tab = left_join(interval.tab, read_tab)
interval.tab$counts = interval.tab$counts %>% tidyr::replace_na(0)
#let's normalize for visualization
interval.tab$probability_scaled = interval.tab$probs_bin %>% scales::rescale(to = c(0,max(interval.tab$counts)))

#plotting
print("plotting...")
ggplot(interval.tab) +
  geom_area(aes(transcriptomic_distance,counts), fill="cornflowerblue", size=0.5) +
  geom_line(aes(transcriptomic_distance,probability_scaled), color="red",size=1) +
  theme_classic() +
  ggtitle("Fragments distribution on transcriptomic space")

# [Writing]
#++++++++++
print("writing...")
ggsave(args$OUTPUT_PDF, scale = 1, device = "pdf",units = "in", width = 11.69, height = 8.27)
fwrite(interval.tab[,c('transcriptomic_distance','counts','probs_bin')],
       file = args$OUTPUT_PROB, sep="\t", col.names = F, nThread=1)

