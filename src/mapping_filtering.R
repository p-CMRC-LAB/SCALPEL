
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))


# Argument Parser
# +++++++++++++++
parser = ArgumentParser(description='Filtering of unmappeds reads and Cleaning')
parser$add_argument('bed', type='character', help='path of bed_gtf file')
parser$add_argument('exons', type='character', help='path of exons file')
parser$add_argument('tr_distance', type='character', help='transcriptomic distance end threshold value')
parser$add_argument('output_path', type='character', help='path of output bed file')
args = parser$parse_args()


#params
TRANSCRIPTOMIC_DISTANCE = as.numeric(args$tr_distance)

#0. Opening
#++++++++++
message("reading ops...")
#reads
reads = fread(
  args$bed,
  col.names = c("seqnames.rd", "start.rd", "end.rd", "strand.rd", "read.id", "frag.id", "splice", "nb.splices"),
  nThread=1)
#exons
exons = fread(
  args$exons, nThread=1)


#encoding read.id & frag.id
message("encoding...")
reads$read.id.encoded = as.numeric(as.factor(reads$read.id))
reads$frag.id.encoded = as.numeric(as.factor(reads$frag.id))

#save
encoded.saves = distinct(reads, read.id, frag.id, read.id.encoded, frag.id.encoded)

#subsetting
reads = reads[,-c("read.id","frag.id")]


#1. Mapping
#++++++++++

#1a. find overlaps
message("overlapping ops...")
hits.gr = findOverlaps(
  query = reads %>% makeGRangesFromDataFrame(start.field = "start.rd", end.field = "end.rd",
                                              strand.field = "strand.rd", seqnames.field = "seqnames.rd"),
  subject = exons %>% makeGRangesFromDataFrame()
)

#2. merging
#++++++++++
message("merging ops...")
merged = cbind(reads[hits.gr@from,], exons[hits.gr@to,-c("seqnames","gene_id","transcript_id")])

#3. filtering ops
#++++++++++++++++
message("filtering ops...")

#- discard fragments associated to intronic mapped reads
message("discard fragments associated to intronic mapped reads")
frag_todel = reads %>%
  dplyr::filter(!(frag.id.encoded %in% unique(merged$frag.id.encoded))) %>%
  distinct(frag.id.encoded)
merged = dplyr::filter(merged, !(frag.id.encoded %in% frag_todel$frag.id.encoded))

#- discard fragments associated to reads overlapping other exons than the 3'end region
message("discard fragments associated to reads overlapping other exons than the 3'end region")
merged$ft.encoded = as.numeric(as.factor(paste0(merged$transcript_name,"***",merged$frag.id.encoded)))
frag_todel = merged %>%
  distinct(transcript_name, frag.id.encoded, read.id.encoded, start.rd, end.rd, start, end, exon_number,ft.encoded) %>%
  group_by(transcript_name) %>%
  filter(exon_number != 1) %>%
  filter(start.rd < start | end.rd > end) %>%
  distinct(ft.encoded) %>%
  data.table()
merged = dplyr::filter(merged, !(ft.encoded %in% frag_todel$ft.encoded))

#- for each fragment, considers the transcript with the highest number of reads
message("for each fragment, considers the transcript with the highest number of reads")
frag_todel = merged %>%
  distinct(gene_name, transcript_name, frag.id.encoded, read.id.encoded, start.rd, end.rd, ft.encoded) %>%
  group_by(frag.id.encoded, transcript_name) %>%
  mutate(n_reads = n()) %>%
  group_by(frag.id.encoded) %>%
  filter(n_reads < max(n_reads)) %>%
  distinct(ft.encoded) %>%
  data.table()
merged = dplyr::filter(merged, !(ft.encoded %in% frag_todel$ft.encoded))
merged$row.id = 1:nrow(merged)

#- check spliced reads coherency
message("check spliced reads coherency")
spliceds = distinct(merged, gene_name, transcript_name, exon_id, exon_number, start, end, 
         frag.id.encoded, read.id.encoded, start.rd, end.rd, strand, nb.splices, ft.encoded, row.id) %>%
  filter(nb.splices>1) %>%
  arrange(read.id.encoded)

unspliceds = distinct(merged, gene_name, transcript_name, exon_id, exon_number, start, end, 
                      frag.id.encoded, read.id.encoded, start.rd, end.rd, strand, nb.splices, ft.encoded, row.id) %>%
  filter(nb.splices==1) %>%
  arrange(read.id.encoded)

#---- bordering criteria
message("bordering criteria")
frag_todel = spliceds %>%
  filter((start.rd != start) & (end.rd != end)) %>%
  distinct(ft.encoded)
spliceds = dplyr::filter(spliceds, !(ft.encoded %in% frag_todel$ft.encoded))
unspliceds = dplyr::filter(unspliceds, !(ft.encoded %in% frag_todel$ft.encoded))

#---- exon consecutivity
message("exon consecutivity")
#a
frag_todel = spliceds %>%
  distinct(transcript_name,read.id.encoded,exon_id,exon_number,nb.splices,ft.encoded) %>%
  group_by(transcript_name,read.id.encoded) %>%
  filter(n()==1) %>%
  data.table() %>%
  distinct(ft.encoded)
spliceds = dplyr::filter(spliceds, !(ft.encoded %in% frag_todel$ft.encoded))
unspliceds = dplyr::filter(unspliceds, !(ft.encoded %in% frag_todel$ft.encoded))

#b
frag_todel = spliceds %>%
  distinct(transcript_name,read.id.encoded,exon_id,exon_number,nb.splices,ft.encoded) %>%
  group_by(transcript_name,read.id.encoded) %>%
  mutate(consecutivity = max(diff(exon_number)), nb_exons = n_distinct(exon_id)) %>%
  filter(consecutivity!=1 | nb_exons != nb.splices) %>%
  data.table() %>%
  distinct(ft.encoded)
spliceds = dplyr::filter(spliceds, !(ft.encoded %in% frag_todel$ft.encoded))
unspliceds = dplyr::filter(unspliceds, !(ft.encoded %in% frag_todel$ft.encoded))

#---- concordance unspliceds/spliceds fragments
message("concordance unspliceds/spliceds fragments")
#(? is there spliceds fragment not associated to the same transcript in unspliceds ?)
tmp = distinct(spliceds, frag.id.encoded, transcript_name, ft.encoded)
frag_todel = dplyr::filter(unspliceds, frag.id.encoded %in% tmp$frag.id.encoded) %>%
  filter(!(ft.encoded %in% tmp$ft.encoded)) %>%
  distinct(ft.encoded)
unspliceds = dplyr::filter(unspliceds, !(ft.encoded %in% frag_todel))

#- merging
reads = merged %>%
  filter(row.id %in% (rbind(spliceds,unspliceds) %>% distinct(row.id))$row.id)

#free memory space
rm(unspliceds)
rm(spliceds)
rm(tmp)
rm(frag_todel)
rm(merged)
gc()


#4. calculates relatives coordinates
#+++++++++++++++++++++++++++++++++++
message("calculate relative coords...")

#POSITIVE strand
tab_pos = reads %>% filter(strand == "+")
tab_pos = tab_pos %>% mutate(rel_start_rd = end_rel + (end - start.rd), rel_end_rd = end_rel + (end-end.rd))
pos_targets = distinct(tab_pos, transcript_name, frag.id.encoded, read.id.encoded, start.rd, end.rd, rel_start_rd, rel_end_rd ,exon_id, strand) %>%
  group_by(transcript_name, frag.id.encoded) %>%
  mutate(rel_start_fg = min(rel_start_rd), rel_end_fg = max(rel_end_rd)) %>%
  data.table() %>%
  distinct(exon_id, frag.id.encoded, read.id.encoded, rel_start_fg, rel_end_fg)
#merging
tab_pos = left_join(tab_pos, pos_targets)
tab_pos$dist_END = tab_pos$rel_end_fg

#NEGATIVE strand
tab_neg = reads %>% filter(strand == "-")
tab_neg = tab_neg %>% mutate(rel_start_rd = end_rel - (end-start.rd), rel_end_rd = end_rel - (end-end.rd))
neg_targets = distinct(tab_neg, transcript_name, frag.id.encoded, read.id.encoded, start.rd, end.rd, rel_start_rd, rel_end_rd ,exon_id, strand) %>%
  group_by(transcript_name,frag.id.encoded) %>%
  mutate(rel_start_fg = min(rel_start_rd), rel_end_fg = max(rel_end_rd)) %>%
  data.table() %>%
  distinct(transcript_name, exon_id, frag.id.encoded, read.id.encoded, rel_start_fg, rel_end_fg)
#merging
tab_neg = left_join(tab_neg, neg_targets)
tab_neg$dist_END = tab_neg$rel_start_fg

#biding tables
reads = data.table(rbind(tab_pos, tab_neg))

#filtering based on distance
targets = (reads %>% filter(dist_END>TRANSCRIPTOMIC_DISTANCE))$ft.encoded
reads = reads %>% filter(!ft.encoded %in% targets)

#column_selection
reads = reads %>%
  dplyr::select(seqnames.rd,start.rd,end.rd,start,end,tr_length,start_rel,end_rel,rel_start_rd,rel_end_rd,
                dist_END,strand,read.id.encoded,frag.id.encoded,splice,nb.splices,gene_name,transcript_name,
                exon_id,exon_number,collapsed_trs,bulk_weights,rel_start_fg,rel_end_fg,ft.encoded)

#free memory space
rm(tab_pos)
rm(pos_targets)
rm(tab_neg)
rm(neg_targets)
rm(targets)
gc()


#5. Writing
#++++++++++
message("writing...")
#filtered table
fwrite(reads, file = args$output_path, sep="\t", row.names = F, nThread=1)
#encoding tables
fwrite(encoded.saves, file = "encoding.txt", sep="\t", row.names = F)


