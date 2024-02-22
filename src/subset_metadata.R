
library(argparse)
library(dplyr)
library(data.table)


#Argument Parser
#===============
parser = ArgumentParser(description='Subset reads info columns')
parser$add_argument('bedf', type="character", help='path of bed file')
parser$add_argument('output_path', type="character", help='path of output bed file')
args = parser$parse_args()

#open file
print("open files...")
bed = fread(args$bedf, col.names= c("seqnames","start.rd","end.rd","strand.rd","read.id","frag.id"), nThread = 1)
# colnames(bed) = c("chr","Start","End","read_id","Strand","bc","umi")

#discard duplicates
bed = bed %>% distinct(seqnames,start.rd,end.rd,strand.rd,frag.id, .keep_all = T)

#add + 1 to the start coord to make it match with the 1 based GTF format
bed$start.rd = bed$start.rd + 1

#check for the scifi notation...
bed$start.rd = format(bed$start.rd, scientific=FALSE)
bed$end.rd = format(bed$end.rd, scientific=FALSE)

#get the max number of splices for each read_id
print("get max reads...")
rsplit = stringr::str_split_fixed(bed$read.id, "/",n=2)
bed$read.id = rsplit[,1]
bed$splice = rsplit[,2]
bed = bed %>% mutate(splice = ifelse(splice=="",0,splice))
bed = bed %>% group_by(read.id) %>% mutate(nb_splices = max(splice)) %>% ungroup()

#arrange reads
bed = bed %>% arrange(seqnames,start.rd,end.rd,strand.rd,read.id,frag.id)

#write
print("writing..")
fwrite(bed, file=args$output_path,sep="\t", row.names=F, col.names=F)





