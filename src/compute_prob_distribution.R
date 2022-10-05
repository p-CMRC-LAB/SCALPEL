

#Libraries
# library(MASS, lib.loc = "/CEPH/shared/R/4.0.2/packages")
# library(survival, lib.loc = "/CEPH/shared/R/4.0.2/packages")
# library(evd)
library(fitdistrplus)
library(dplyr)
library(data.table)
library(ggplot2)
# library(ggthemes, lib.loc = "/CEPH/shared/R/4.0.2/packages")
library(argparse)
library(scales)
library(distributionsrd, lib.loc = '~/CEPH/RPACKAGES/')
# library("kdensity")

#Variables
QUANTILE_THRESHOLD = "60%"
BINS = 25



#Argunent Parser
#**************
parser = ArgumentParser(description='Probability of unique transcript')
parser$add_argument('PATH_OF_UNIQUE_TR', metavar='Tfip', type="character", help='path of unique isoform all reads')
parser$add_argument('OUTPUT_PROB', metavar='Tfop', type="character", help='path of output unique isoform intervals probabilities')
parser$add_argument('OUTPUT_PDF', metavar='PDFfop', type="character", help='path of output pdf pictures')
args = parser$parse_args()


#Files Opening (1)
#*************
#get list of of unique read file
# files = list.files('~/CEPH/scalpel_rev/fd/fd_files/reads/fragments/', pattern = "*.fragment_filtered_unique", full.names = T)
files = list.files(args$PATH_OF_UNIQUE_TR, pattern = "*.fragment_filtered_unique", full.names = T)
reads = parallel::mclapply(files, function(x) fread(x), mc.preschedule = T, mc.cores = 5)
reads = do.call(rbind, reads)
reads

#filter out genes bringing variabilities
genes_counts = table(reads$gene_name)
genes_quants = genes_counts %>% quantile(seq(0,1,0.01))
reads = reads %>% filter(gene_name %in% names(genes_counts[genes_counts<genes_quants[[QUANTILE_THRESHOLD]]]))
reads

#get distinct readid/3end
dtab = reads %>% distinct(read_id,dist_END) %>% arrange(dist_END)
dtab = table(dtab$dist_END) %>% data.table()
dtab$V1 = as.numeric(dtab$V1)
colnames(dtab) = c('transcriptomic_distance','counts')

#create a table covering the whole transcriptomic space
MIN = min(reads$dist_END)
MAX = max(reads$dist_END)
dtab_cov = data.table(data.frame(transcriptomic_distance = seq(MIN,MAX))); dtab_cov
#merge
merged = left_join(dtab_cov, dtab)
merged$counts = merged$counts %>% tidyr::replace_na(0); dtab
#eliminate duplicates
merged = merged %>% distinct(.keep_all=TRUE)

#plot
ggplot(merged, aes(transcriptomic_distance,counts)) +
  geom_point() +
  geom_line() +
  geom_smooth(span=0.30) +
  theme_bw()

#breaks
breaks = unique(c(-Inf,sort(unique(c(seq(0,MIN,-BINS),seq(0,MAX,BINS)))),MAX)); breaks
merged.2 = do.call(rbind, lapply(1:(length(breaks)-1), function(x){
  A = breaks[x]
  B = breaks[x+1]
  ends = merged$transcriptomic_distance[(merged$transcriptomic_distance>A) & (merged$transcriptomic_distance<=B)]
  ends_counts = merged$counts[merged$transcriptomic_distance %in% ends]
  return(data.frame(transcriptomic_distance = ends,bin_counts = sum(ends_counts)))
}))
merged.2 = merged.2 %>% distinct(transcriptomic_distance, .keep_all=TRUE)

#final_merging
merged.final = left_join(merged.2, merged) %>% data.table(); merged.final
# merged.final$densityk = density(merged.final$counts, n = length(merged.final$counts))$y %>% rescale(to=c(0,1))

#final scaling
merged.final$bin_counts_pb = merged.final$bin_counts / sum(merged.final$bin_counts)
merged.final$probability_normalized = merged.final$bin_counts_pb %>% rescale(to=c(0,1))
merged.final$counts_normalized = merged.final$counts %>% rescale(to=c(0,1))
merged.final = merged.final[,c('transcriptomic_distance','counts','bin_counts_pb','counts_normalized','probability_normalized','bin_counts')]
colnames(merged.final) = c('transcriptomic_distance','counts','probs_bin','counts_normalized','probability_normalized','bin_counts')
merged.final$probs_bin = merged.final$probs_bin * 100
merged.final


merged.loess = loess(probs_bin~transcriptomic_distance, data = merged.final[merged.final$transcriptomic_distance>=0], span = 0.5)
merged.final$loess = 0
merged.final$loess[merged.final$transcriptomic_distance>=0] = merged.loess$fitted
merged.final


fd = fitdist((merged.final$bin_counts/1), 'frechet', start = list(shape=1,scale=1))
merged.final$fd = dfrechet(0:(nrow(merged.final)-1), shape = fd$estimate[1], scale = fd$estimate[2]) %>% rescale(to=c(0,1))


#plot
ggplot(merged.final) +
  geom_line(aes(transcriptomic_distance,probability_normalized, col='empiric_distribution')) +
  geom_line(aes(transcriptomic_distance,loess, col='loess')) +
  # geom_smooth(aes(transcriptomic_distance, probability_normalized), span=0.2) +
  # geom_line(aes(transcriptomic_distance, densityk, col='density'), stat = 'identity') +
  # geom_line(aes(transcriptomic_distance, fd, col='frechet'), stat = 'identity') +
  # geom_point(aes(transcriptomic_distance,counts_normalized)) +
  geom_line(aes(transcriptomic_distance,counts_normalized, col='read_counts')) +
  theme_bw(base_size = 20)

#writing
merged.final$loess = merged.final$loess * 100
merged.final$probs_bin = merged.final$loess
ggsave(args$OUTPUT_PDF, scale = 2)
fwrite(merged.final[,c('transcriptomic_distance','counts','probs_bin','counts_normalized','probability_normalized')], file = args$OUTPUT_PROB, sep="\t")






# df_tab = df_tab[,c('dist_END','counts','probs_bin','counts_normalized','probability_normalized')]
# descdist(merged$counts, discrete = FALSE)
#
#
# #let's orient it on a [0_Inf axis]
# df$dist_END_adapted = df$dist_END + (abs(min(df$dist_END)) + 1)
# descdist(df$dist_END, discrete = FALSE)
# df2 = table(df$dist_END) %>% data.table(); df2
# df2$loess = (loess(N~V1, data=df2, span=0.3))$fitted
# ggplot(df2) +
#   geom_line(aes(as.numeric(V1),N)) +
#   geom_line(aes(as.numeric(V1),loess), col='red',size=1) +
#   geom_point(aes(as.numeric(V1),N))
# colnames(df2) = c('dist_END','counts','probs_bin')
#
#
#
#
#
#
#
#
#
#
# # # betaf = fitdist(rescale(df$dist_END, to=c(0,1)), 'beta')
# # # dbeta(min(df$dist_END_adapted):max(df$dist_END_adapted), shape1 = 1.13, shape2 = 1.20)
# # #
# # # #hist
# # # bins = 10
# # # breaks = c(seq(0,min(df$dist_END)-bins,-bins),seq(0,max(df$dist_END)+bins,bins)) %>% sort %>% unique()
# # # res = hist(df$dist_END, breaks = breaks)
# # # # t1 = data.frame(breaks = res$breaks[2:length(res$breaks)], counts=res$counts) %>% data.table()
# # # # fd = fitdist(t1$counts, 'frechet',start=list(shape=1,scale=5), method = 'mge'); fd
# # # # t1$fd = dfrechet(1:nrow(t1),shape=fd$estimate[1], scale = fd$estimate[2])
# # # # t1
# # # # t1$cn = t1$counts %>% rescale(to=c(0,1))
# # # # t1$fd = t1$fd %>% rescale(to=c(0,1))
# # # #
# # # # ggplot(t1) +
# # # #   geom_line(aes(breaks,cn)) +
# # # #   geom_point(aes(breaks,fd))
# # #
# # # #let's fit a frechet distribution
# # # frechet_fit = fitdist((df$dist_END_adapted)/1,'frechet',start=list(shape=1,scale=5), method = 'mge'); frechet_fit
# # # nbinom_fit = fitdist((df$dist_END_adapted)/1,'nbinom', method = 'mge'); nbinom_fit
# # # poisson_fit = fitdist((df$dist_END_adapted)/1,'pois'); poisson_fit
#
# df_counts = data.table(df %>% group_by(dist_END) %>% summarize(read_id=n()))
# df_counts$dist_END = as.character(df_counts$dist_END)
#
# df_tab = data.table()
# df_tab$dist_END = min(df$dist_END):max(df$dist_END)
# df_tab$dist_END = as.character(df_tab$dist_END)
#
# # # df_tab$frechet_dstb = dfrechet(min(df$dist_END_adapted):max(df$dist_END_adapted), shape = frechet_fit$estimate[1], scale = frechet_fit$estimate[2])
# # # df_tab$nbinom_dstb = dnbinom(min(df$dist_END_adapted):max(df$dist_END_adapted), size = nbinom_fit$estimate[1], mu = nbinom_fit$estimate[2])
# # # df_tab$npoisson_dstb = dpois(min(df$dist_END_adapted):max(df$dist_END_adapted), lambda = poisson_fit$estimate[1])
# #
#
# df_tab = left_join(df_tab, df2, on='dist_END'); df_tab
# df_tab$counts[is.na(df_tab$counts)] = 0
# df_tab$probs_bin = (loess(counts~dist_END, data = df_tab))$fitted
# df_tab$dist_END = as.numeric(df_tab$dist_END)
# # colnames(df_tab) = c('dist_END','probs_bin','counts')
#
#
# #scaling
# df_tab$counts_normalized = df_tab$counts %>% rescale(to=c(0,1))
# # df_tab$probability_normalized = df_tab$probs_bin %>% rescale(to=c(0,1))
# df_tab$probability_normalized = (df_tab$probs_bin / sum(df_tab$probs_bin)) %>% rescale(to=c(0,1))
#
# pt = ggplot(df_tab) +
#   # geom_line(aes(dist_END, loess)) +
#   geom_line(aes(dist_END,counts_normalized, color='obs')) +
#   geom_line(aes(dist_END,probability_normalized, color='nbinom_distribution'), size=2) +
#   ggtitle("reads distribution on transcriptomic space") +
#   theme_bw(base_size = 15); pt
#  #writing
# ggsave(args$OUTPUT_PDF, device = 'jpeg',scale = 0.7)
# df_tab = df_tab[,c('dist_END','counts','probs_bin','counts_normalized','probability_normalized')]
# fwrite(df_tab, file = args$OUTPUT_PROB, sep="\t")
