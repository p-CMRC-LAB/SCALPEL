
suppressMessages(suppressWarnings(library(argparse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(GenomicRanges)))


# Argument Parser
# +++++++++++++++
parser = ArgumentParser(description='Filtering of unmappeds reads and Cleaning')
parser$add_argument('bed', type='character', help='path of bed_gtf file')
parser$add_argument('exons', type='character', help='path of exons file')
parser$add_argument('output_path', type='character', help='path of output bed file')
args = parser$parse_args()


#0. Opening
#----------
message("reading ops...")
#reads
reads = fread(
  args$bed,
  col.names = c("seqnames.rd","start.rd","end.rd","strand.rd","read.id","frag.id"),
  nThread=1)
#exons
gtf = fread(
  args$exons, nThread=1)


#BED file 1-based conversion
reads$start.rd = reads$start.rd + 1


#Processing of reads
#------------------
message("Processing of reads...")
reads = reads %>%
  dplyr::filter(start.rd>min(gtf$start)-1e5 & end.rd<max(gtf$end)+1e5) %>%
  #discard pcr replicates
  distinct(seqnames.rd,start.rd,end.rd,strand.rd,frag.id, .keep_all = T)

#Mapping reads into the genome
#-----------------------------
message("Mapping reads into the genome...")
hits = GenomicRanges::findOverlaps(
  GenomicRanges::makeGRangesFromDataFrame(reads, seqnames.field = "seqnames.rd", start.field = "start.rd",
                                          end.field = "end.rd", strand.field = "strand.rd"),
  GenomicRanges::makeGRangesFromDataFrame(gtf, seqnames.field = "seqnames", start.field = "start",
                                          end.field = "end", strand.field = "strand")
)
mapped = cbind(reads[hits@from,], gtf[hits@to,]) %>%
  dplyr::filter(strand.rd==strand) %>% 
  dplyr::select(!seqnames) %>%
  tidyr::separate(col="read.id", into=c("readID","splice.pos"), sep="/", remove = F) %>%
  data.table() %>%
  dplyr::mutate(splice.pos=as.numeric(as.character(splice.pos)))


#Filtering operations
#--------------------
message("Filtering operations...")

#a. discard fragments associated to intergenics/intronics reads
message("a. discard fragments associated to intergenics/intronics reads...")
reads = mapped %>%
  filter(!frag.id %in% (reads %>% filter(!(read.id %in% mapped$read.id)))$frag.id)

rm(mapped)
rm(hits)
gc()

#b. discard fragments associated to reads overlapping transcripts coordinates
message("b. discard fragments associated to reads overlapping transcripts coordinates...")
trs.todel = (reads %>%
               dplyr::mutate(check=ifelse(exon_number==1 & strand=="+" & start.rd<start,"wrong_overlap","correct"),
                      check=ifelse(exon_number==1 & strand=="-" & end.rd>end,"wrong_overlap",check),
                      check=ifelse(exon_number!=1 & (start.rd<start | end.rd>end), "wrong_overlap",check)) %>%
               dplyr::filter(check=="wrong_overlap") %>%
               distinct(frag.id, transcript_name) %>%
               dplyr::reframe(ftrs = paste0(frag.id,"_",transcript_name)))$ftrs
reads = reads %>%
  dplyr::mutate(ftrs = paste0(frag.id,"_",transcript_name)) %>%
  dplyr::filter(!(ftrs %in% trs.todel))

#c. For each fragment, check if the mapped transcripts intersect all the associated read ids
message("c. For each fragment, check if the mapped transcripts intersect all the associated read ids...")
trs.todel = (reads %>%
  group_by(frag.id) %>%
  mutate(nb.readIDs_tot = n_distinct(readID)) %>%
  group_by(frag.id,transcript_name) %>%
  dplyr::reframe(nb.readIDs_counted = n_distinct(readID), nb.readIDs_tot,ftrs) %>%
  distinct() %>%
  dplyr::filter(nb.readIDs_tot!=nb.readIDs_tot) %>%
  data.table())$ftrs
reads = reads %>%
  dplyr::filter(!(ftrs %in% trs.todel))


#d. checks spliced reads coherency
message("d. checks spliced reads coherency...")
reads = reads %>%
  group_by(readID) %>%
  mutate(nb.splices = n_distinct(splice.pos), splice.pos = as.numeric(as.character(splice.pos)),
         splice.pos = ifelse(is.na(splice.pos), 0, splice.pos)) %>%
  data.table()

spliceds = reads %>% dplyr::filter(nb.splices>1) %>% arrange(readID)
unspliceds = reads %>% dplyr::filter(nb.splices==1) %>% filter(splice.pos<=1)

#---- bordering criteria
message("bordering criteria...")
trs.todel = (spliceds %>%
  dplyr::filter((start.rd != start) & (end.rd != end)) %>%
  distinct(ftrs))$ftrs
spliceds = spliceds %>% dplyr::filter(!(ftrs %in% trs.todel))
unspliceds = unspliceds %>% dplyr::filter(!(ftrs %in% trs.todel))

#---- Number of splices on a transcript 
message("nb splices coherency...")
#a
trs.todel = (spliceds %>%
  group_by(transcript_name,readID) %>%
  dplyr::filter(n_distinct(read.id)!=nb.splices) %>% distinct(ftrs))$ftrs
spliceds = dplyr::filter(spliceds, !(ftrs %in% trs.todel))
unspliceds = dplyr::filter(unspliceds, !(ftrs %in% trs.todel))

#---- exon consecutivity coherency
message("exon consecutivity...")
trs.todel = (spliceds %>%
  group_by(transcript_name,readID) %>%
  mutate(consecutivity = abs(max(diff(exon_number))), nb_exons = n_distinct(exon_number)) %>%
  dplyr::filter(consecutivity!=1 | nb_exons != nb.splices) %>%
  distinct(ftrs))$ftrs
spliceds = dplyr::filter(spliceds, !(ftrs %in% trs.todel))
unspliceds = dplyr::filter(unspliceds, !(ftrs %in% trs.todel))


#---- concordance unspliced/spliced fragments
message("concordance unspliceds/spliceds fragments")
#(? is there spliced reads of a fragment not associated to the same transcript in spliced/unspliced ?)
tmp = distinct(spliceds, frag.id, transcript_name, ftrs)
trs.todel = dplyr::filter(unspliceds, frag.id %in% tmp$frag.id) %>%
  dplyr::filter(!(ftrs %in% tmp$ftrs)) %>%
  distinct(ftrs)
unspliceds = dplyr::filter(unspliceds, !(ftrs %in% trs.todel))

#--- Merging unspliceds and spliceds
reads = rbind(spliceds, unspliceds) %>%
  dplyr::select(!c(read.id,splice.pos,ftrs)) %>%
  distinct()

#free memory space
rm(unspliceds)
rm(spliceds)
rm(tmp)
rm(trs.todel)
gc()

print(reads)

#Calculate relative coordinates
#------------------------------
message("Calculate relative coordinates...")
reads = reads %>%
  dplyr::mutate(start.rdR = ifelse(strand=="+",endR-(end.rd-start),NA),
         end.rdR = ifelse(strand=="+",endR-(start.rd-start),NA),
         start.rdR = ifelse(strand=="-",endR-(end-start.rd),start.rdR),
         end.rdR = ifelse(strand=="-",endR-(end-end.rd),end.rdR))

#Writing
#-------
message("writing...")
reads %>% 
  dplyr::select(!c(collapsed,nb.splices,collapsed)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(transcript_name, frag.id) %>%
  dplyr::mutate(fg.start = min(start.rd), fg.end=max(end.rd)) %>%
  ungroup() %>%
  saveRDS(file = args$output_path)
  
