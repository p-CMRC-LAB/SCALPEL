

#!/usr/bin/env Rscript


#********************************
#Function to compute EM algorithm
#********************************


#libraries
#*********
library(data.table)
library(dplyr)
#library(reshape, lib.loc = "~/CEPH/R_PACKAGES/")
library(reshape)
library(argparse)

#settings
#********
options(future.globals.maxSize = 1500 * 1024^2)
ITERATIONS_MAX = 100


#Argunent Parser
#**************
parser = ArgumentParser(description='Comput EM algorithm... ')
parser$add_argument('PATH_OF_CELL_FILE', metavar='PRFfip', type="character",
                    help='path of input cell probs')
parser$add_argument('PATH_OF_PRED_FILE', metavar='PPFfop', type="character",
                    help='path of output predicted cells')
args = parser$parse_args()

#function posterior
posterior = function(probabilities, theta_values){
    numerateur = probabilities * theta_values
    denominateur = (probabilities * theta_values)
    return(numerateur / sum(denominateur))
}


#Opening
#*******
# get probability file
cat("reading...\n")
prob_file = fread(args$PATH_OF_CELL_FILE)
# prob_file = fread('~/CEPH/fresh_vs_ad_experiment/ad/scalpel/RESULTS/READS/CELLS/AAAACAATAGGA.cell')

#convert probabilities to numeric
prob_file$probs_bin = as.numeric(prob_file$probs_bin)
weights = data.frame(prob_file %>% distinct(gene_name,transcript_name,salmon_rlp))


#list of uniques GENES
GENES = prob_file$gene_name %>% unique()
#Final table
FINAL = prob_file %>% dplyr::select(bc,gene_name,gene_id,transcript_name,transcript_id) %>% distinct()
FINAL$probs_bin = 0

cat("looping\n")
GENES = sort(as.character(GENES))
for(GENE in GENES){
    print(GENE)
    #read by genes
    A = prob_file %>% filter(gene_name == GENE)
    Wgene = weights %>% filter(gene_name == GENE)
    rownames(Wgene) = Wgene$transcript_name
    # print(A)
    # print(Wgene)

    #case 1 fragment 1 Isoform
    if(nrow(A) == 1){
        FINAL[which(FINAL$gene_name == GENE), "probs_bin"] = 1
    }else{
        A = A %>% cast(umi~transcript_name, fill = 0);
        #case several fragments but 1 isoform
        if(ncol(A) == 2){
            FINAL[which(FINAL$gene_name == GENE), "probs_bin"] = 1
            next
        }
        #a - get initial isoforms probabilities
        A = A[,2:ncol(A)];
        # print(A)

        #multilply fragments probabililies by the transcripts probability quantified by salmon
        cltemp = colnames(A)
        A = (t(t(A) * Wgene[colnames(A),'salmon_rlp']))
        colnames(A) = cltemp
        # print(A)

        estimated_abundances = rep(1/ncol(A),ncol(A))
        previous_abundances = estimated_abundances

        for(i in 1:ITERATIONS_MAX){
            #calculate posterior
            outA = apply(A,1, function(x) {posterior(x, estimated_abundances)});
            estimated_abundances = apply(outA, 1, mean)
            if( !(FALSE %in% ((abs(previous_abundances - estimated_abundances) < 0.001))) ){
                break
            }
            previous_abundances = estimated_abundances
        }
        FINAL[which(FINAL$gene_name == GENE), "probs_bin"] = round(estimated_abundances,3)
    }
}

cat("writing...\n")
fwrite(FINAL, file = args$PATH_OF_PRED_FILE, quote =F, sep = "\t", row.names = F, col.names = T)
