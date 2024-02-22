
library(argparse)
suppressWarnings(library(data.table))
suppressWarnings(library(dplyr))
suppressWarnings(library(stringr))
suppressWarnings(library(GenomicRanges))


# Argument Parser
# -----------------
parser = ArgumentParser(description='Filtering of unmappeds reads and Cleaning')
parser$add_argument('bed', type='character', help='path of bed_gtf file')
parser$add_argument('exons', type='character', help='path of exons file')
parser$add_argument('tr_distance', type='character', help='transcriptomic distance end threshold value')
parser$add_argument('output_path', type='character', help='path of output bed file')
args = parser$parse_args()

#params
TRANSCRIPTOMIC_DISTANCE = as.numeric(args$tr_distance)

#Reads
print("opening...")
bed.input = fread(args$bed, col.names = c("seqnames", "start.rd", "end.rd", "strand.rd", "read.id", "frag.id", "splice", "nb_splices"), nThread = 1)
bed.gr = bed.input %>% GenomicRanges::makeGRangesFromDataFrame(start.field = "start.rd", end.field = "end.rd", strand.field = "strand.rd", seqnames.field = "seqnames", keep.extra.columns = T)
#exons
exons.input = fread(args$exons, nThread = 1)
exons.gr =  exons.input %>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T)



#Mapping
#=======
print("read/exons mapping...")
#mapping
hits.gr = findOverlaps(bed.gr, exons.gr) %>% data.frame()
bed = data.frame(bed.gr)
colnames(bed) = c("seqnames.rd","start.rd","end.rd","width.rd","strand.rd","read.id","frag.id","splice","nb.splices")
bed = bed[hits.gr$queryHits,] %>% data.table()
exons = data.frame(exons.gr)[hits.gr$subjectHits,] %>% data.table()

#c.binding
reads = cbind(bed,exons)

#discard all the frags associated to outmapped reads (1)
targets = bed.input %>% filter(!read.id %in% reads$read.id) %>% distinct(frag.id)
reads = reads %>% dplyr::filter(!frag.id %in% targets$frag.id)
reads = reads %>% filter(strand == strand.rd)

#renames columns
reads = reads %>% dplyr::select(seqnames,start.rd,end.rd,start,end,width.rd,tr_length,start_rel,end_rel,strand,read.id,frag.id,splice,nb.splices,gene_id,gene_name,transcript_id,
                                transcript_name,exon_id,exon_number,collapsed_trs,bulk_weights)
colnames(reads) = c("seqnames","start.rd","end.rd","start","end","width.rd","tr_length","start_rel","end_rel","strand","read.id","frag.id","splice","nb.splices","gene_id","gene_name",
                    "transcript_id","transcript_name","exon_id","exon_number","collapsed_trs","bulk_weights")

print("splicing/filtering...")
#get for each fragment the max number of reads associated on a transcript (2)
reads$fidt = paste0(reads$frag.id, "_", reads$transcript_name)
targets = reads %>%
    distinct(frag.id,transcript_name,read.id) %>%
    group_by(frag.id, transcript_name) %>%
    summarise(read.nb=n()) %>%
    group_by(frag.id) %>%
    filter(read.nb!=max(read.nb)) %>% ungroup()
targets$fidt = paste0(targets$frag.id, "_", targets$transcript_name)
reads = reads %>% filter(!fidt %in% targets$fidt)


#Discard all the fragments associated to reads mapped out of transcripts (2bis)
reads$fidt = paste0(reads$frag.id, "_", reads$transcript_name)
targets = reads %>% 
    dplyr::select(start.rd,end.rd,start,end,strand,exon_number,transcript_name,fidt) %>%
    mutate(check_notlast = ifelse(exon_number!=1 & (start.rd<start | end.rd>end), TRUE,FALSE),
           check_lastpos = ifelse((exon_number==1 & strand=="+" & start.rd<start), TRUE,FALSE),
           check_lastneg = ifelse((exon_number==1 & strand=="-" & end.rd>end), TRUE,FALSE)) %>%
    filter(check_notlast==T | check_lastpos==T | check_lastneg==T)
reads = reads %>% filter(!fidt %in% targets$fidt)


#Get for each fragment the most mapped transcript (3)
targets = reads %>%
    distinct(frag.id,transcript_name,exon_id) %>%
    group_by(frag.id, transcript_name) %>%
    summarise(exon.nb=n()) %>%
    group_by(frag.id) %>%
    filter(exon.nb!=max(exon.nb)) %>% ungroup()
targets$fidt = paste0(targets$frag.id, "_", targets$transcript_name)
reads = reads %>% filter(!fidt %in% targets$fidt)

#get spliceds reads
spliceds = reads %>% filter(nb.splices>1)
unspliceds = reads %>% filter(nb.splices==1)

#filter based on bordering (1)
targets = (spliceds %>% filter(start.rd!=start & end.rd!=end))$fidt
spliceds = spliceds %>% filter(!fidt %in% targets)
unspliceds = unspliceds %>% filter(!fidt %in% targets)

#check exon consecutivity (2)
targets = spliceds %>%
    distinct(read.id, transcript_name, exon_id, exon_number, nb.splices, fidt) %>%
    group_by(transcript_name,read.id) %>%
    mutate(check1 = length(unique(diff(exon_number))), check2 = length(unique(exon_id))) %>%
    filter(check1!=1 | check2!=nb.splices) %>%
    ungroup()
spliceds = spliceds %>% filter(!fidt %in% targets$fidt)
unspliceds = unspliceds %>% filter(!fidt %in% targets$fidt)

#check concordance frag spliceds/unspliceds (4)
check.tab = spliceds %>% dplyr::select(frag.id, transcript_name, fidt) %>% distinct()
# extract the reads for which the frag_ID is present in the unspliced & spliced table, but not associated with the transcript it should be in the unspliced table based on the spliced table
targets = (unspliceds %>% filter(frag.id %in% check.tab$frag.id & (!fidt %in% check.tab$fidt)))$fidt 
unspliceds = unspliceds %>% filter(!fidt %in% targets)

#merge tables
reads = rbind(spliceds,unspliceds)
rm(unspliceds)
rm(spliceds)


#calculates relatives coordinates
#================================
print("calculate relative coords...")
#pos strand
tab_pos = reads %>% filter(strand == "+")
tab_pos = tab_pos %>% mutate(rel_start_rd = end_rel + (end - start.rd), rel_end_rd = end_rel + (end-end.rd))
tab_pos$dist_END = tab_pos$rel_end_rd
#neg_strand
tab_neg = reads %>% filter(strand == "-")
tab_neg = tab_neg %>% mutate(rel_start_rd = end_rel - (end-start.rd), rel_end_rd = end_rel - (end-end.rd))
tab_neg$dist_END = tab_neg$rel_start_rd
reads = rbind(tab_pos, tab_neg)

#filtering based on distance
targets = (reads %>% filter(dist_END>TRANSCRIPTOMIC_DISTANCE))$fidt
reads = reads %>% filter(!fidt %in% targets)

#select columns
reads = reads %>% dplyr::select(seqnames,start.rd,end.rd,start,end,width.rd,tr_length,start_rel,end_rel,rel_start_rd,rel_end_rd,dist_END,strand,read.id,frag.id,splice,nb.splices,gene_id,gene_name,transcript_id,transcript_name,exon_id,exon_number,collapsed_trs,bulk_weights,fidt)


#Writing
#=======
fwrite(reads, file=args$output_path, sep="\t", row.names=F)