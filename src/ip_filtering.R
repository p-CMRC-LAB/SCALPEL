

suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))


# Argument Parser
# +++++++++++++++
parser = ArgumentParser(description='Filtering of internal priming reads')
parser$add_argument('reads', type='character', help='path of bed_gtf file')
parser$add_argument('ip_file', type="character", help='path of internal priming ref file')
parser$add_argument('threshold_dist', type="character", help='distance of ip position from isoform end')
parser$add_argument('output_read_ids', type="character", help='path of output file for read ids mapped')
parser$add_argument('output_bed_unique', type="character", help='path of output file for read mapped in unique transcript gene')
parser$add_argument('output_ip', type="character", help='path of output ips')
parser$add_argument('output_path', type='character', help='path of output reads file')
args = parser$parse_args()

#params
THRESHOLD_IP_DIST = as.numeric(args$threshold_dist)

#0. Opening
#-----------
message("Opening...")
#reads
reads = readRDS(args$reads)
#ipdb
ipdb = fread(args$ip_file,col.names = c("seqnames.ip","start.ip","end.ip","c4","c5","strand.ip"), nThread=1)

if(nrow(ipdb)==0){

    #return void ipdb
    message("write ip sites mapped...")
    ipdb %>% fwrite(file=args$output_ip, sep="\t",row.names=F,col.names=F, nThread=2)

    #return reads
    #filtering
    reads = reads %>% 
     distinct(seqnames.rd,start.rd,end.rd,strand.rd,start.rdR,frag.id,start,end,gene_name,transcript_name,bulk_TPMperc,readID)

    #writing readID
    distinct(reads, readID) %>% fwrite(file=args$output_read_ids, sep="\t")
    reads = reads %>% dplyr::select(!readID)


    #4. Filter reads associated to genes with unique transcripts
    unique_genes = (reads %>% distinct(gene_name,transcript_name) %>% 
     group_by(gene_name) %>%
     summarise(count_tr = n()) %>%
     dplyr::filter(count_tr == 1))$gene_name
    reads_unique = reads %>% dplyr::filter(gene_name %in% unique_genes) %>% distinct()


    #5. Writing reads table
    message("writing...")
    reads = reads %>% dplyr::rename(dist_END="start.rdR")
    saveRDS(reads, file=args$output_path)
    #fwrite(reads, file=args$output_path, sep="\t", row.names=F, col.names=T, nThread=1)
    fwrite(reads_unique, file=args$output_bed_unique, sep="\t", row.names=F, col.names=F,  nThread=1)


}else {

    ipdb = ipdb %>% dplyr::select(!c(c4,c5))

    #Mapping
    #-------
    message("mapping...")
    annots = distinct(reads, seqnames.rd,start,end,startR,endR,strand,gene_name,transcript_name) ; annots
    hits = GenomicRanges::findOverlaps(
      GenomicRanges::makeGRangesFromDataFrame(annots, seqnames.field = "seqnames.rd", keep.extra.columns = T),
      GenomicRanges::makeGRangesFromDataFrame(ipdb, seqnames.field = "seqnames.ip", start.field = "start.ip",
                                              end.field = "end.ip", strand.field = "strand.ip")
    )
    message("internal priming processing...")
    ipdb = cbind(annots[hits@from,], ipdb[hits@to,]) %>%
      dplyr::filter(start.ip>start & end.ip<end) %>%
      #calculation of internal priming relative coordinates
      mutate(start.ipR = ifelse(strand=="-",endR-(end-start.ip),NA),
             end.ipR = ifelse(strand=="-",endR-(end-end.ip),NA),
             start.ipR = ifelse(strand=="+",endR-(end.ip-start),start.ipR),
             end.ipR = ifelse(strand=="+",endR-(start.ip-start),end.ipR)) %>%
      dplyr::filter(start.ipR>=THRESHOLD_IP_DIST) %>%
      #Get the closest IP out of threshold distance for each transcript
      group_by(transcript_name) %>%
      dplyr::filter(start.ipR == min(start.ipR)) %>%
      distinct(transcript_name,gene_name,start.ip,end.ip) %>%
      data.table(); ipdb

    #3. Write internal priming annotation sites
    message("write ip sites mapped...")
    ipdb %>% fwrite(file=args$output_ip, sep="\t",row.names=F,col.names=F, nThread=2)


    #Merging
    message("ip sites filtering...")
    trs_todel = (reads %>%
      left_join(ipdb) %>%
      #discard all fragment upstream of an IP position
      mutate(check=ifelse((!is.na(start.ip)) & strand=="+" & fg.end<start.ip,T,F),
             check=ifelse((!is.na(start.ip)) & strand=="-" & fg.start>end.ip,T,F),
             ftrs = paste0(transcript_name,"_",frag.id)) %>%
      dplyr::filter(check==T) %>%
      distinct(ftrs))$ftrs

    #filtering
    reads = reads %>%
      mutate(ftrs=paste0(transcript_name,"_",frag.id)) %>%
      dplyr::filter(!ftrs %in% trs_todel) %>%
      distinct(seqnames.rd,start.rd,end.rd,strand.rd,start.rdR,frag.id,start,end,gene_name,transcript_name,bulk_TPMperc,readID)

    #writing readID
    distinct(reads, readID) %>% fwrite(file=args$output_read_ids, sep="\t")
    reads = reads %>% dplyr::select(!readID)


    #4. Filter reads associated to genes with unique transcripts
    unique_genes = (reads %>% distinct(gene_name,transcript_name) %>%
                      group_by(gene_name) %>%
                      summarise(count_tr = n()) %>%
                      dplyr::filter(count_tr == 1))$gene_name
    reads_unique = reads %>% dplyr::filter(gene_name %in% unique_genes) %>% distinct()


    #5. Writing reads table
    message("writing...")
    reads = reads %>% dplyr::rename(dist_END="start.rdR")
    saveRDS(reads, file=args$output_path)
    fwrite(reads_unique, file=args$output_bed_unique, sep="\t", row.names=F, col.names=F,  nThread=1)

}
