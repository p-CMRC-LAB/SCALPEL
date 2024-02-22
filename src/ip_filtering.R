

library(argparse)
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(GenomicRanges))


# Argument Parser
# -----------------
parser = ArgumentParser(description='Filtering of internal priming reads')
parser$add_argument('bed', type='character', help='path of bed_gtf file')
parser$add_argument('ip_file', type="character", help='path of internal priming ref file')
parser$add_argument('threshold_dist', type="character", help='distance of ip position from isoform end')
parser$add_argument('output_read_ids', type="character", help='path of output file for read ids mapped')
parser$add_argument('output_bed_unique', type="character", help='path of output file for read mapped in unique transcript gene')
parser$add_argument('output_ip', type="character", help='path of output ips')
parser$add_argument('output_path', type='character', help='path of output bed file')
args = parser$parse_args()


#params
THRESHOLD_IP_DIST = as.numeric(args$threshold_dist)

#reads
reads = fread(args$bed, nThread =1)
#ipdb
ipdb = fread(args$ip_file, col.names = c("seqnames","start.ip","end.ip","c4","c5","strand.ip"), nThread=1)


#Filter Internal priming sites
#==============================
print("internal priming filtering...")

if(length(ipdb)!=0){
    #annotate ip positions with isoforms infos
    ipdb.gr = ipdb[,c("seqnames","start.ip","end.ip","strand.ip")] %>% GenomicRanges::makeGRangesFromDataFrame(start.field = "start.ip", end.field = "end.ip", strand.field = "strand.ip")
    #merge with exon annotation
    exons.gr.mapped = reads %>% dplyr::distinct(seqnames,start,end,strand,transcript_name,gene_name,exon_id,exon_number,start_rel,end_rel,tr_length) %>%
        GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)
    hits = findOverlaps(ipdb.gr,exons.gr.mapped) %>% data.frame()
    ipdb = cbind(ipdb = data.frame(ipdb.gr)[hits$queryHits,], data.frame(exons.gr.mapped)[hits$subjectHits,]) %>% data.table()
    ipdb = ipdb %>% filter(ipdb.strand == strand)
    ipdb = ipdb[,-c("strand","seqnames")]
    
    #filter ip_positions
    #calculate rel coords ipdb
    ipdb = ipdb %>% filter(ipdb.start>=start & ipdb.end<=end)
    #pos
    ip_pos = ipdb %>% filter(ipdb.strand=="+")
    ip_pos = ip_pos %>% mutate(rel_start_ip = end_rel + (end-ipdb.start), rel_end_ip = end_rel + (end - ipdb.end))
    ip_pos$dist_END_ip = ip_pos$rel_end_ip
    #neg
    ip_neg = ipdb %>% filter(ipdb.strand=="-")
    ip_neg = ip_neg %>% mutate(rel_start_ip = end_rel - (end - ipdb.start), rel_end_ip = end_rel - (end-ipdb.end))
    ip_neg$dist_END_ip = ip_neg$rel_start_ip
    ipdb = rbind(ip_pos,ip_neg) %>% distinct()
    #get the closest ip out of threshold for each transcrip
    ipdb = ipdb %>% filter(dist_END_ip > THRESHOLD_IP_DIST)
    ipdb = ipdb %>% filter(transcript_name %in% reads$transcript_name)
    
    #filtering
    targets = ipdb %>% group_by(transcript_name) %>% filter(dist_END_ip == min(dist_END_ip)) %>% ungroup()
    targets = left_join(reads, targets[,c("rel_start_ip","rel_end_ip","transcript_name","gene_name","dist_END_ip")], multiple="all") %>%
        na.omit() %>%
        filter((strand=="+" & rel_end_rd>dist_END_ip) | (strand=="-" & rel_start_rd>dist_END_ip))
    reads$fidt = paste0(reads$frag.id, "_", reads$transcript_name)
    reads = reads %>% filter(!fidt %in% targets$fidt)
    reads = reads[,-c("fidt")]
}

print("writing...")
#get read_IDS to filter
reads %>% distinct(read.id) %>% fwrite(args$output_read_ids,sep="\t",row.names = F,col.names = F)

#write internal priming annotation sites
ipdb %>% fwrite(file=args$output_ip, sep="\t",row.names=F,col.names=F, nThread=2)

# Filter reads associated to genes with unique transcripts (5)
unique_genes = (reads %>% distinct(gene_name,transcript_name) %>% 
    group_by(gene_name) %>%
    summarise(count_tr = n()) %>%
    filter(count_tr == 1))$gene_name
reads_unique = reads %>% filter(gene_name %in% unique_genes) %>% distinct()

#write file
fwrite(reads, file=args$output_path, sep="\t", row.names=F, col.names=T, nThread=2)
fwrite(reads_unique, file=args$output_bed_unique, sep="\t", row.names=F, col.names=F,  nThread=2)

