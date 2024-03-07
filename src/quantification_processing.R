


#author: Franz AKE
#July 2022

#libs
library(argparse)
library(data.table)
library(dplyr)

#Argument parser
parser = ArgumentParser(description='Processing of salmon quantification_file')
parser$add_argument('salmon_quant', type="character", help='path of salmon *.sf file')
parser$add_argument('genome', type="character", help='path of gtf file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()



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


qf = dplyr::distinct(A, gene_name, transcript_id, transcript_name, Length, NumReads)

#writing
fwrite(qf, sep="\t", row.names=F, file =args$output_path, nThread=1)