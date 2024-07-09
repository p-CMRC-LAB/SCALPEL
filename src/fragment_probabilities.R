
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))

#Argument Parser
#+++++++++++++++
parser = ArgumentParser(description='Probability of fragments')
parser$add_argument('reads_path', type="character", help='path of  all regads')
parser$add_argument('probs_path', type="character", help='path of probs')
<<<<<<< HEAD
=======
parser$add_argument('encoding', type="character", help='path of encoding gfile')
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
parser$add_argument('output_path', type="character", help='output')
args = parser$parse_args()


#0. Opening
#++++++++++
message("Opening...")
<<<<<<< HEAD
reads = readRDS(args$reads_path)
probs = fread(args$probs_path, col.names = c("dist_END","counts","probs_bin"))
=======
reads = fread(args$reads_path, nThread =1)
probs = fread(args$probs_path, col.names = c("dist_END","counts","probs_bin"), nThread =1)
#encode.tab
encode.tab = fread(args$encoding)


>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

#1. Joining
#++++++++++
reads = left_join(reads, probs)
reads$probs_bin = replace_na(reads$probs_bin, 0)
reads$counts = replace_na(reads$counts, 0)
<<<<<<< HEAD
#discard reads with null probability
reads = dplyr::filter(reads, probs_bin!=0)

#2. Probability weighting and column selection
#+++++++++++++++++++++++++++++++++++++++++++++
reads = reads %>% dplyr::distinct(start.rd,end.rd,frag.id,gene_name,transcript_name,bulk_TPMperc,probs_bin)

#splitting
reads = reads %>% tidyr::separate(frag.id, into = c("bc","umi"), sep="::")
=======


#2. Probability weighting and column selection
#+++++++++++++++++++++++++++++++++++++++++++++
reads = reads %>% dplyr::distinct(start.rd,end.rd,read.id.encoded,frag.id.encoded,gene_name,transcript_name,
                                  bulk_weights,probs_bin)
reads$probs_weighted = reads$probs_bin

#decode frag.id
reads = left_join(
  reads,
  encode.tab %>% 
    filter(frag.id.encoded %in% reads$frag.id.encoded) %>%
    distinct(frag.id.encoded,frag.id))

#splitting
reads = reads %>% tidyr::separate(frag.id, into = c("bc","umi"), sep="::") %>% data.table()
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

#3. Grouping
#+++++++++++
reads = reads %>% group_by(bc,umi,gene_name,transcript_name) %>% 
<<<<<<< HEAD
  mutate(frag_probs = prod(probs_bin)) %>% #intersection of reads probabilities ! are there independant events ? yes
  dplyr::distinct(bc,gene_name,transcript_name,umi,frag_probs,bulk_TPMperc) %>%
  mutate(frag_probs_weighted = frag_probs * bulk_TPMperc) %>%
  arrange(bc,gene_name,transcript_name,umi) %>%
  dplyr::distinct(bc,gene_name,transcript_name,umi,frag_probs_weighted) %>%
  ungroup()
=======
  mutate(frag_probs_weighted = prod(probs_weighted)) %>% #intersection of reads probabilities ! are there independant events ? yes
  dplyr::select(bc,gene_name,transcript_name,umi,frag_probs_weighted,bulk_weights) %>%
  mutate(frag_probs_weighted = frag_probs_weighted * bulk_weights) %>%
  arrange(bc,gene_name,transcript_name,umi) %>%
  dplyr::select(bc,gene_name,transcript_name,umi,frag_probs_weighted) %>%
  distinct() %>% 
  data.table()
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

#delete sequencing tags
reads$bc = str_replace(reads$bc, pattern = "XC:Z:","")
reads$bc = str_replace(reads$bc, pattern = "CB:Z:","")
reads$umi = str_replace(reads$umi, pattern = "XM:Z:","")
reads$umi = str_replace(reads$umi, pattern = "UB:Z:","")

#4. Writing
#++++++++++
<<<<<<< HEAD
fwrite(reads, file=args$output_path, sep="\t", col.names=F, row.names=F)
=======
fwrite(reads, file=args$output_path, sep="\t", col.names=F, row.names=F, nThread=1)
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47

