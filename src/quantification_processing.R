
#libs
library(argparse)
library(dplyr)

#Argument parser
parser = ArgumentParser(description='Processing of salmon quantification_file')
parser$add_argument('genome', type="character", help='path of gtf file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()


#Opening
#GTF
gtf = rtracklayer::import(args$genome) %>%
    data.frame() %>%
    dplyr::filter(type=="exon")

#SAMPLES QUANTIFICATION
lapply(list.files(path = ".", pattern = "*.sf", full.names = T), function(x){
    #reading bulk counts
    curr = data.table::fread(x, nThread=1)
    curr$sample = x
    return(curr)
}) %>% data.table::rbindlist() %>%
    #get MeanTPM & discard null transcrips in bulk
    group_by(Name) %>%
    summarise(meanTPM = mean(TPM)) %>%
    dplyr::filter(meanTPM!=0) %>%
    rename(transcript_id="Name") %>%
    #include gene metadata
    left_join(distinct(gtf, gene_name, transcript_name, transcript_id)) %>%
    #calculate Transcripts bulk % in gene
    group_by(gene_name) %>%
    dplyr::summarise(transcript_name, bulk_TPMperc = meanTPM / sum(meanTPM)) %>%
    arrange(gene_name,transcript_name) %>%
    data.table::fwrite(file=args$output_path, nThread=1, row.names=F, sep="\t")

#SPLIT GTF by chromosomes
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
    x = distinct(x, seqnames,start,end,width,strand,gene_name,transcript_name,exon_number)
    data.table::fwrite(x, file = paste0(x$seqnames[1],".gtf"), sep="\t", nThread=1)
})
