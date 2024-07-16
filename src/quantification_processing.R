
#libs
suppressWarnings(library(argparse, quietly=T, verbose=F))
suppressWarnings(library(dplyr, quietly=T, verbose=F))

#Argument parser
parser = ArgumentParser(description='Processing of salmon quantification_file')
parser$add_argument('genome', type="character", help='path of gtf file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()


#Opening

#GTF
message("GTF opening...")
gtf = rtracklayer::import(args$genome) %>%
    data.frame() %>%
    dplyr::filter(type=="exon")

message("GTF check the expected attributes to be present...")
if( FALSE %in% c(c("gene_id","transcript_id") %in% colnames(gtf)) ){ stop("Error! The attributes 'gene_id' or 'transcript_id' are expected in your GTF file") }

if( ("gene_name" %in% colnames(gtf)) & ("transcript_name" %in% colnames(gtf)) ) {
    gtf = gtf %>%
        dplyr::mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name), transcript_name = ifelse(is.na(transcript_name), transcript_id, transcript_name)) %>%
        dplyr::distinct(seqnames, start, end, width, strand, gene_id, gene_name, transcript_id, transcript_name)
} else {
    gtf = gtf %>%
        dplyr::mutate(gene_name = gene_id, transcript_name = transcript_id) %>%
        dplyr::distinct(seqnames, start, end, width, strand, gene_id, gene_name, transcript_id, transcript_name)
}

#check1
if(nrow(gtf)==0){ stop("Error! GTF file contains no 'exon' entries ! Please provide valid GTF file") }

#SAMPLES QUANTIFICATION
message("Processing bulk quantification sample files...")
lapply(list.files(path = ".", pattern = "*.sf", full.names = T), function(x){
    #reading bulk counts
    curr = data.table::fread(x, nThread=1)
    curr$sample = x
    return(curr)
}) %>% data.table::rbindlist() %>%
    #get MeanTPM & discard null transcrips in bulk
    dplyr::group_by(Name) %>%
    dplyr::summarise(meanTPM = mean(TPM)) %>%
    dplyr::filter(meanTPM!=0) %>%
    dplyr::rename(transcript_id="Name") %>%
    #include gene metadata
    left_join(distinct(gtf, gene_name, transcript_name, transcript_id)) %>%
    #calculate Transcripts bulk % in gene
    dplyr::group_by(gene_name) %>%
    dplyr::reframe(transcript_name, bulk_TPMperc = meanTPM / sum(meanTPM)) %>%
    arrange(gene_name,transcript_name) -> bulk_quants

#writing merge quants
bulk_quants %>% stats::na.omit() %>% data.table::fwrite(file=args$output_path, nThread=1, row.names=F, sep="\t")

#filtering of GTF
gtf = dplyr::filter(gtf, gene_name %in% bulk_quants$gene_name) %>% stats::na.omit()

#check2
if(nrow(gtf)==0){ stop("Error! No gene in the GTF have been bulk quantified by Salmon quant! Please provide valid GTF file") }

#SPLIT GTF by chromosomes
message("GTF splitting...")
lapply(split(gtf, gtf$seqnames), function(x){
    #in case of GTF file from Ensembl

    # a. rename columns
    lookup = c(gene_type = "gene_biotype", transcript_type = "transcript_biotype")
    if("gene_biotype" %in% colnames(x)){
        rename(x, all_of(lookup))
    }

    # b. check gene version
    if("gene_version" %in% colnames(x)){
        x$gene_id = paste0(x$gene_id, ".", x$gene_version)
        x$transcript_id = paste0(x$transcript_id, ".", x$transcript_version)
    }

    #writing
    x = distinct(x, seqnames,start,end,width,strand,gene_name,transcript_name)
    if(nrow(x)!=0){ 
        data.table::fwrite(x, file = paste0(x$seqnames[1],".gtf"), sep="\t", nThread=1)
    }
}) -> writeds
