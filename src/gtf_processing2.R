

library(argparse)
library(dplyr)
library(data.table)
library(stringr)

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
DIST_TR = as.numeric(args$DIST_TR)
DIST_TR_EXON = as.numeric(args$DIST_TR_EXON)


#Function to calculate relative coordinates
#==========================================
rel_coord = function(vec, strand){
    start_rel = c()
    end_rel = c()
    out = list()
    for(i in 1:length(vec)){
        if(i==1){
            if(strand=="-"){
                start_rel = append(start_rel,0)
                end_rel = append(end_rel, vec[i]-1)
            }else{
                start_rel = append(start_rel,vec[i]-1)
                end_rel = append(end_rel, 0)
            }
        }else{
            if(strand=="-"){
                start_rel = append(start_rel, vec[i-1])
                end_rel = append(end_rel, (vec[i-1] + vec[i]))
            }else{
                start_rel = append(start_rel, (vec[i-1] + vec[i]))
                end_rel = append(end_rel, vec[i-1])
            }
            vec[i] = (vec[i-1] + vec[i]) + 1
        }
    }
    #return
    out$starts = start_rel
    out$ends = end_rel
    return(out)
}

#Function for calculation of relative exon coordinates
#=====================================================
get_rel_coords = function(width, strand){
    rel_coords = rel_coord(width, strand[1])
    return(paste(rel_coords$starts, rel_coords$ends, sep = "_"))
}


#Function for collapsing isoforms
#================================
collapse_isoforms = function(tab, dist_tr, dist_cov){
    # Function for collapsing isofrms with similar ends
    # =================================================
    strand_in = tab$strand[1]
    
    if(nrow(tab)==1){return(NULL)}
    #Let's loop into the exons (1)
    for(ind in sort(unique(tab$exon_number))){
        
        curr_tab = tab %>% filter(exon_number == ind)
        
        if(ind == 1){
            #if last exon, look for isoforms with similar ends...
            if(strand_in == "+"){
                curr_tab$cluster_ends = MsCoreUtils::group(curr_tab$end, dist_tr)
                curr_tab$cluster_ends = paste0(curr_tab$cluster_ends, "_", curr_tab$start_rel)
                curr_tab = curr_tab %>% group_by(cluster_ends) %>%
                    mutate(tr2 = paste0(transcript_name, collapse = "_")) %>%
                    filter(transcript_name != tr2)
                tab = tab %>% filter(transcript_name %in% curr_tab$transcript_name)
            }else{
                curr_tab$cluster_ends = MsCoreUtils::group(curr_tab$start, dist_tr)
                curr_tab$cluster_ends = paste0(curr_tab$cluster_ends, "_", curr_tab$end_rel)
                curr_tab = curr_tab %>% group_by(cluster_ends) %>%
                    mutate(tr2 = paste0(transcript_name, collapse = "_")) %>%
                    filter(transcript_name != tr2)
                tab = tab %>% filter(transcript_name %in% curr_tab$transcript_name)
            }
            
            #break
            if(ind == max(tab$exon_number)){break}
            if(nrow(curr_tab)==0){return(NULL)}
            
        }else{
            if(strand_in == "+"){
                curr_tab = curr_tab %>% group_by(end) %>%
                    mutate(start = ifelse(start_rel==dist_cov, end - (start_rel - max(end_rel)),start)) %>% ungroup()
            }else{
                curr_tab = curr_tab %>% group_by(start) %>%
                    mutate(end = ifelse(end_rel==dist_cov, start + (end_rel - max(start_rel)),end)) %>% ungroup()
            }
            curr_tab = curr_tab %>% group_by(start,end) %>%
                mutate(tr2 = paste0(transcript_name, collapse="_")) %>%
                filter(transcript_name != tr2)
            tab = tab %>% filter(transcript_name %in% curr_tab$transcript_name)
            
            #break
            if(nrow(tab)== 0){break}
            if(nrow(curr_tab)==0){return(NULL)}
        }
    }
    
    #table of isoforms with similar ends
    return(unique(tab$transcript_name))
}



# 1) File opening
#================
print("File opening...")
gtf = fread(GTF_PATH)
qf = fread(QF_PATH)


# 2) Process annotation
#======================
print("process annotation...")

# a. rename columns
lookup = c(gene_type = "gene_biotype", transcript_type = "transcript_biotype")
if("gene_biotype" %in% colnames(gtf)){
    rename(gtf, all_of(lookup))
}

# b. check gene version
if("gene_version" %in% colnames(gtf)){
    gtf$gene_id = paste0(gtf$gene_id, ".", gtf$gene_version)
    gtf$transcript_id = paste0(gtf$transcript_id, ".", gtf$transcript_version)
}

# c. Integrate salmon Bulk quantification of isoforms (check why Salmon does not provide all the transcript even with a null value)
#discard all the transcript with a null quantification at Bulk level
print(gtf)
gtf = left_join(gtf, dplyr::select(qf, transcript_id, NumReads)) %>% 
    filter(!is.na(NumReads)) %>%
    filter(NumReads != 0)

# c. exon selection & relative coordinates calculation
print("relative exon coordinates calculation...")
print(gtf)
gtf = gtf %>%
    filter(type=="exon") %>%
    dplyr::select(seqnames,start,end,strand,width,type,gene_id,gene_name,gene_type,
                  transcript_id,transcript_name,transcript_type,exon_number,exon_id,NumReads) %>%
    mutate(exon_number = as.numeric(exon_number)) %>%
    arrange(seqnames,gene_name,transcript_name,desc(exon_number)) %>%
    group_by(transcript_name) %>%
    mutate(rel_coords = get_rel_coords(width,strand), tr_length = sum(width))
rel_coords = str_split_fixed(gtf$rel_coords, pattern = "\\_", n=2)
gtf$start_rel = as.numeric(rel_coords[,1])
gtf$end_rel = as.numeric(rel_coords[,2])
gtf = dplyr::select(gtf, !(rel_coords)) %>% ungroup()
gtf = gtf %>% arrange(transcript_name,exon_number) %>% group_by(transcript_name) %>%
    mutate(exon_number = rev(exon_number))

#TRUNCATING isoforms based on distance threshold
gtf_neg = gtf %>% filter(strand == "+" & end_rel < DIST_TR) %>%
    group_by(transcript_name) %>%
    mutate(start = ifelse(start_rel>DIST_TR, (end-(DIST_TR-end_rel)), start), start_rel = ifelse(start_rel>DIST_TR,DIST_TR,start_rel)) %>%
    data.table()
gtf_pos = gtf %>% filter(strand == "-" & start_rel < DIST_TR) %>%
    group_by(transcript_name) %>%
    mutate(end = ifelse(end_rel>DIST_TR, (start+(DIST_TR-start_rel)), end), end_rel = ifelse(end_rel>DIST_TR,DIST_TR,end_rel)) %>%
    data.table()
gtf = rbind(gtf_neg, gtf_pos)


# e. Get similar isoforms in the distance threshold range
print("get similar isoforms...")
RES = lapply(split(gtf, gtf$gene_name), function(x) collapse_isoforms(x, dist_tr = DIST_TR_EXON, dist_cov = DIST_TR))
sim_gtf = gtf %>% filter(transcript_name %in% unlist(RES)) %>%
    distinct(transcript_name, gene_name) %>%
    group_by(gene_name) %>%
    mutate(tr_clusters = paste0(transcript_name, collapse="_")) %>%
    data.table()

# Filtering...
gtf = left_join(gtf, sim_gtf)
gtf = gtf %>%
    mutate(tr_clusters = ifelse(is.na(tr_clusters), transcript_name, tr_clusters)) %>%
    group_by(tr_clusters) %>%
    filter(NumReads == max(NumReads)) %>%
    data.table()

# f. Calculate average weigth TPM for each isoforms / gene
gtf = gtf %>% group_by(gene_name) %>%
 mutate(salmon_rlp = NumReads / sum(unique(NumReads))) %>%
 mutate(salmon_rlp = round(salmon_rlp,2)) %>%
 dplyr::filter(salmon_rlp != 0) %>% 
 data.table()
gtf = subset(gtf, select = -c(NumReads))


# g. ordering
gtf = gtf %>% arrange(gene_name,transcript_name,exon_number)

# h. Get unique genes associated to chr
gtf_unique = gtf %>% dplyr::select(gene_name,transcript_name) %>% 
distinct() %>% 
group_by(gene_name) %>%
 mutate(counts = n()) %>%
  filter(counts == 1) %>%
   ungroup()
gtf_unique = gtf_unique[,c("gene_name","counts")]
fwrite(gtf_unique, file=args$OUTPUT_UNIQUE, sep="\t")

# i. select columns and rename
# gtf$exon_number = gtf$exon_ind
gtf = gtf %>%
    dplyr::select(seqnames,start,end,strand,start_rel,end_rel,tr_length,gene_id,
                  gene_name,transcript_id,transcript_name,exon_id,exon_number,
                  tr_clusters,salmon_rlp)
colnames(gtf) = c("seqnames","start","end","strand","start_rel","end_rel","tr_length",
                  "gene_id","gene_name","transcript_id","transcript_name","exon_id",
                  "exon_number","collapsed_trs","bulk_weights")

# j. ordering
gtf = gtf %>% arrange(seqnames,start,end,strand)

# h. write processed GTF into a file
print("writing...")
fwrite(gtf, file = args$OUTPUT, sep="\t")
