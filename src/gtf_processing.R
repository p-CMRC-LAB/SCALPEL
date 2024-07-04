

library(argparse)
library(dplyr)
library(data.table)

#ArgumentParser
#==============
parser = ArgumentParser(description='Process GTF file')
parser$add_argument('GTF_PATH', type="character", help='path of gtf file')
parser$add_argument('QF_PATH', type="character", help='path of salmon file')
parser$add_argument('DIST_TR', type="character", help='transcriptomic distance range threshold')
parser$add_argument('DIST_TR_EXON', type="character", help='transcripts inter distance threshold')
parser$add_argument('OUTPUT', type="character", help='path of output gtf file')
parser$add_argument('OUTPUT_UNIQUE', type="character", help='path of output gtf file unique isoform')
args = parser$parse_args()


GTF_PATH = args$GTF_PATH
QF_PATH = args$QF_PATH
DT_THRESHOLD = as.numeric(args$DIST_TR)
DT_EX_THRESHOLD = as.numeric(args$DIST_TR_EXON)


#Functions
#----------

groupExon_ends = function(x, distance_exon){
  #Function for grouping exon with similar ends
  #=============================================
  if (is.unsorted(x)) {
    idx <- order(x)
    x <- x[idx]
  } else idx <- integer()
  res <- cumsum(c(1L, diff(x) > distance_exon))
  res[idx] <- res
  res
}

collapseIsoforms = function(x, distance_transcriptome, distance_ex){
  #function for collapsing Isoforms
  #================================

  #get isoforms with similar ends
  strand = x$strand[1]
  x.tmp = x %>%
    dplyr::filter(exon_number==1)
  x.tmp = x.tmp %>%
    dplyr::mutate(clusters = ifelse(strand=="-",groupExon_ends(start, distance_ex),groupExon_ends(end, distance_ex)))
  simIsoforms = group_by(x.tmp, clusters) %>% dplyr::filter(n()>1)
  simIsoforms = lapply(split(simIsoforms, simIsoforms$clusters), function(x) x$transcript_name)
  x = left_join(x, distinct(x.tmp, transcript_name, clusters))

  #filter input tab
  collapsed = data.frame(transcript_name=character(0),collapsed=character(0))
  #loop on similar isoforms
  for(trs in simIsoforms){
    curr = dplyr::filter(x, transcript_name %in% trs)
    for(ind in unique(curr$exon_number)){
      current = curr %>% dplyr::filter(exon_number==ind)
      #in last exons
      if(ind==1){
        if(strand=="-"){
          current = mutate(current, end= ifelse(endR==distance_transcriptome-1,max(end),end))
          current.2 = group_by(current, clusters, end) %>%
            mutate(collapsed = paste(transcript_name, collapse = "_"), nCollapsed = n_distinct(transcript_name)) %>%
            dplyr::filter(nCollapsed!=1) %>% ungroup()
        }else{
          current = mutate(current, start= ifelse(endR==distance_transcriptome-1,min(start),start))
          current.2 = group_by(current, clusters, start) %>%
            mutate(collapsed = paste(transcript_name, collapse = "_"), nCollapsed = n_distinct(transcript_name)) %>%
            dplyr::filter(nCollapsed!=1) %>% ungroup()
        }
        #breaking rule
        if(nrow(current.2)==0){break}
        if(ind == max(curr$exon_number)){
          collapsed = rbind(collapsed, distinct(current.2, transcript_name, collapsed))
        }
      }else{
        if(strand=="-"){
          current.2 = group_by(current,clusters,start) %>%
            mutate(end= ifelse(n()>1 & endR==distance_transcriptome-1,max(end),end)) %>% ungroup()
        }else{
          current.2 = group_by(current,clusters,end) %>%
            mutate(start= ifelse(n()>1 & endR==distance_transcriptome-1,min(start),start)) %>% ungroup()
        }
        current.2 = group_by(current,start,end) %>%
          mutate(collapsed=NA, collapsed = paste(transcript_name, collapse="_"), nCollapsed = n_distinct(transcript_name)) %>%
          dplyr::filter(nCollapsed!=1) %>% ungroup()
        #breaking rule
        if(nrow(current.2)==0){break}
        if(ind == max(curr$exon_number)){
          collapsed = rbind(collapsed, distinct(current.2, transcript_name, collapsed))
        }
      }
    }
  }
  return(stats::na.omit(collapsed))
}



# 1) File opening
#----------------
print("File opening...")
gtf = fread(GTF_PATH)
qf = fread(QF_PATH, col.names=c("gene_name","transcript_name", "bulk_TPMperc")) %>% dplyr::filter(transcript_name %in% gtf$transcript_name)

if(nrow(qf)==0){

    warning(paste("Any bulk quantification found for the file ", GTF_PATH))

}else{

    #Formatting table (A)
    #----------------
    print("Formatting table...")
    gtf = gtf %>%
    distinct(seqnames,start,end,width,strand,gene_name,transcript_name) %>%
    arrange(seqnames,start,desc(end)) %>%
    stats::na.omit() %>%
    group_by(transcript_name) %>%
    dplyr::mutate(exon_number = ifelse(strand=="+",n():1,1:n())) %>%
    arrange(transcript_name) %>%
    data.table()


    #Calculate relative coordinates (B)
    #------------------------------
    print("Calculate relative coordinates...")
    gtf = gtf %>%
    arrange(transcript_name,exon_number) %>%
    group_by(transcript_name) %>%
    group_modify(~{
        starts = .x$start
        ends = .x$end
        startR = rep(NA,length(starts))
        endR = rep(NA,length(ends))
        for(i in .x$exon_number){
        if(i==1){
            startR[i] = 0
            endR[i] = (ends[i] - (starts[i])) - 1
        }else{
            startR[i] = endR[i-1] + 1
            endR[i] = startR[i] + (ends[i] - starts[i])
        }
        .x$startR = startR
        .x$endR = endR
        }
        .x
    }) %>%
    distinct(seqnames,start,end,strand,startR,endR,gene_name,transcript_name,exon_number) %>%
    data.table()


    #Truncate isoforms within transcriptomic distance threshold
    #----------------------------------------------------------
    gtf = gtf %>%
    dplyr::filter(startR < DT_THRESHOLD) %>%
    mutate(start = ifelse(strand == "+" & endR>DT_THRESHOLD,end-(DT_THRESHOLD-startR),start),
            end = ifelse(strand== "-" & endR>DT_THRESHOLD,start+(DT_THRESHOLD-startR),end),
            endR = ifelse(endR>DT_THRESHOLD,DT_THRESHOLD-1,endR)) %>%
    arrange(seqnames,start,desc(end),transcript_name)


    #Collapse similar isoforms
    #-------------------------
    #processing
    collapseds.tab = gtf %>%
    group_by(gene_name) %>%
    group_modify(~{collapseIsoforms(.x, DT_THRESHOLD, DT_EX_THRESHOLD)}) %>%
    ungroup()
    gtf = left_join(gtf, collapseds.tab) %>%
    mutate(collapsed = ifelse(is.na(collapsed),"none",collapsed))

    #writing
    distinct(gtf, transcript_name, collapsed) %>%
    dplyr::filter(collapsed!="none") %>%
    distinct() %>%
    fwrite(file = paste0(gtf$seqnames[1],"_collapsed_isoforms.txt"), sep="\t", col.names = F, row.names = F)

    #collapsing
    gtf = left_join(gtf, qf) %>%
    dplyr::filter(!is.na(bulk_TPMperc)) %>%
    group_by(collapsed) %>%
    mutate(check=ifelse(collapsed!="none",max(bulk_TPMperc),bulk_TPMperc),
            end=ifelse(strand=="+" & collapsed!="none" & exon_number==1,max(end),end),
            start=ifelse(strand=="-" & collapsed!="none" & exon_number==1,min(start),start)) %>%
    dplyr::filter(bulk_TPMperc==check) %>%
    group_by(collapsed) %>%
    mutate(endR=ifelse(collapsed!="none" & exon_number==max(exon_number), (startR + (end-start)) - 1, endR)) %>%
    ungroup() %>%
    distinct(seqnames,start,end,startR,endR,strand,gene_name,transcript_name,exon_number,bulk_TPMperc,collapsed)

    #6) Writing
    #all exons entries
    fwrite(gtf, file = args$OUTPUT, sep="\t", nThread=1)

    #all genes with unique isoforms
    gtf_unique = gtf %>% dplyr::select(gene_name,transcript_name) %>%
    distinct() %>%
    group_by(gene_name) %>%
    mutate(counts = n()) %>%
    dplyr::filter(counts == 1) %>%
    ungroup()
    gtf_unique = gtf_unique[,c("gene_name","counts")]
    fwrite(gtf_unique, file=args$OUTPUT_UNIQUE, sep="\t")

}
