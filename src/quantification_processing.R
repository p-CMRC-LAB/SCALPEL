
<<<<<<< HEAD
#libs
library(argparse)
=======


#author: Franz AKE
#July 2022

#libs
library(argparse)
library(data.table)
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
library(dplyr)

#Argument parser
parser = ArgumentParser(description='Processing of salmon quantification_file')
<<<<<<< HEAD
=======
parser$add_argument('salmon_quant', type="character", help='path of salmon *.sf file')
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
parser$add_argument('genome', type="character", help='path of gtf file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()


<<<<<<< HEAD
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
=======

#opening

#genome
genome = rtracklayer::import(args$genome)
genome = data.table(data.frame((genome)))
genome = genome %>% filter(type=="exon")

#salmon quant
qf = fread(args$salmon_quant, sep ="\t", col.names = c("Name", "Length", "EffectiveLength","TPM", "NumReads", "Sample"))

#expand gene_name
if(stringr::str_detect(qf$Name[1], pattern = "|")){
    qf = qf %>% tidyr::separate(Name, sep = "\\|", into = as.character(c("transcript_id", 1:8))) %>%
        dplyr::select(transcript_id, Length, TPM, Sample) %>%
        data.table()
}else{
    qf = dplyr::select(qf, Name, Length, TPM, Sample)
    colnames(qf)[1] = "transcript_id"
}

#filter null transcripts
qf = distinct(qf)
qf = left_join(qf, distinct(genome, transcript_id, transcript_name, gene_name))
qf = qf %>% group_by(transcript_id) %>%
    mutate(NumReads = round(mean(TPM),10)) %>%
    # filter(NumReads != 0) %>%
    data.table()

# A = qf %>% group_by(gene_name,Sample) %>%
#     mutate(frac_TPM_Sample = TPM / sum(TPM)) %>%
#     data.table()


qf = dplyr::distinct(qf, gene_name, transcript_id, transcript_name, Length, NumReads)

#writing
fwrite(qf, sep="\t", row.names=F, file =args$output_path, nThread=1)
>>>>>>> 0c3fe80469feb05951ba0efe3f2fb8270c284a47
