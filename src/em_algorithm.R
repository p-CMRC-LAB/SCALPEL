


#Libraries
library(dplyr)
library(data.table)
library(argparse)
library(stringr)
library(tidyr)


#Argument Parser
#===============
parser = ArgumentParser(description='Probability of unique transcript')
parser$add_argument('cell', type="character", help='path of cell file')
parser$add_argument('output_path', type="character", help='output path')
args = parser$parse_args()

em_algorithm = function(tab){
    #"""
    # Function to perform EM algorithm...
    # ====================================
    #""""
    #get initial estimated abundances
    tr_size = length(unique(tab$transcript_name))
    estimated_abund = rep(1/tr_size, tr_size)
    res = list()
    res[[1]] = estimated_abund
    intab = dcast(tab, umi ~ transcript_name, fill = 0, value.var = "frag_prob_weighted")
    intab = intab[,2:ncol(intab)]
    for(i in 2:50){
        if(ncol(intab)==1){break}
        #Compute M-step
        current = t(apply(intab,1,function(x) (x*estimated_abund)/sum(x*estimated_abund))) %>% data.table() %>% mutate_all(~replace_na(.,0))
        #Compute E-step
        estimated_abund = apply(current,2,function(x) mean(x))
        res[[i]] = estimated_abund
        #stop criteria for the EM is we encounter a variation of estimated abundance < 0.001
        stop_crit = abs(res[[i]] - res[[i-1]]) %>% max()
        if(stop_crit<0.001){break}
    }
    names(estimated_abund) = colnames(intab)
    db = data.table()
    db$transcript_name = names(estimated_abund)
    db$rel_abund = estimated_abund
    db$gene_name = tab$gene_name[1]
    return(db)
}

#file opening
cell = fread(args$cell, col.names = c("bc","gene_name","gene_id","transcript_name","transcript_id","umi","frag_prob_weighted"), nThread =1) %>% distinct()

#Perform EM algorithm
print("EM algorithm...")
genetabs = split(cell, cell$gene_name)
res = Reduce(rbind, lapply(genetabs, function(x) em_algorithm(x)))

#join with cell information
cell = left_join(cell,res)
cell = cell %>% dplyr::select(bc,gene_name,gene_id,transcript_name,transcript_id,rel_abund)

#writing
fwrite(cell,file = args$output_path, sep="\t", col.names = F, row.names = F, nThread=1)

