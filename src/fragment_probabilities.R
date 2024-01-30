



#Libraries
library(dplyr)
library(data.table)
library(argparse)
library(stringr)
library(tidyr)



#Argunent Parser
#**************
parser = ArgumentParser(description='Probability of fragments')
parser$add_argument('reads_path', type="character", help='path of  all regads')
parser$add_argument('probs_path', type="character", help='path of probs')
parser$add_argument('output_path', type="character", help='output')
args = parser$parse_args()

#file opening
reads = fread(args$reads_path)
probs = fread(args$probs_path, col.names = c("dist_END","counts","probs_bin"))

#joining
reads = left_join(reads, probs, by = c("dist_END"))
reads$probs_bin = replace_na(reads$probs_bin, 0)
reads$counts = replace_na(reads$counts, 0)

#probability weighting and column selection
reads = reads %>% dplyr::distinct(read.id,frag.id,gene_id,gene_name,transcript_id,transcript_name,bulk_weights,probs_bin)
reads$probs_weighted = reads$probs_bin
reads = reads %>% tidyr::separate(frag.id, into = c("bc","umi"), sep="::") %>% data.table()
print(reads)

#grouping
reads = reads %>% group_by(bc,umi,gene_name,transcript_name) %>% 
    mutate(frag_probs_weighted = prod(probs_weighted)) %>% #intersection of reads probabilities ! are there independant events ?
    dplyr::select(bc,gene_name,gene_id,transcript_name,transcript_id,umi,frag_probs_weighted,bulk_weights) %>%
    mutate(frag_probs_weighted = frag_probs_weighted * bulk_weights) %>%
    arrange(bc,gene_name,transcript_name,umi) %>%
    dplyr::select(bc,gene_name,gene_id,transcript_name,transcript_id,umi,frag_probs_weighted) %>%
    distinct() %>% 
    data.table()

print(reads)

#delete sequencing tags
reads$bc = str_replace(reads$bc, pattern = "XC:Z:","")
reads$bc = str_replace(reads$bc, pattern = "CB:Z:","")
reads$umi = str_replace(reads$umi, pattern = "XM:Z:","")
reads$umi = str_replace(reads$umi, pattern = "UB:Z:","")

#writing
fwrite(reads, file=args$output_path, sep="\t", col.names=F, row.names=F)