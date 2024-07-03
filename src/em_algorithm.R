


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


MAX_IT=30

em_algorithm = function(tab){
    #"""
    # Function to perform EM algorithm...
    # ====================================
    #""""
    db = data.table()
    db$transcript_name = character(0)
    db$rel_abund = numeric(0)
    db$gene_name = character(0)
    #get initial estimated abundances
    res = list()
    intab = dcast(tab, umi ~ transcript_name, fill = 0, value.var = "frag_prob_weighted")
    # TODO: Remove zero rows and issue a warning!
    intab = intab[,2:ncol(intab)]
    rsums = rowSums(intab)
    if(all(rsums==0)){
      warning(paste0("All Fragment probabilities of gene", tab$gene_name[1], "are 0 ! Gene discarded !!!"))
      return(db)
    }
    if(any(rsums==0)){
	warning("Fragment with 0 probability in all isoforms ignored!!!")
	intab = intab[rsums!=0,]
    }
    #get initial estimated abundances
    tr_size = length(colnames(intab))
    estimated_abund = rep(1/tr_size, tr_size)
    res[[1]] = estimated_abund
    for(i in 2:MAX_IT){
        if(ncol(intab)==1){break}
        #Compute M-step
        current = t(apply(intab,1,function(x) (x*estimated_abund)/sum(x*estimated_abund))) %>% data.table()# %>% mutate_all(~replace_na(.,0))
        #Compute E-step
        estimated_abund = apply(current,2,function(x) mean(x))
        res[[i]] = estimated_abund
        #stop criteria for the EM is we encounter a variation of estimated abundance < 0.01
        stop_crit = abs(res[[i]] - res[[i-1]]) %>% max()
        if(stop_crit<0.01){break}
    }
    names(estimated_abund) = colnames(intab)
    db = data.table()
    db$transcript_name = names(estimated_abund)
    db$rel_abund = estimated_abund
    db$gene_name = tab$gene_name[1]
    return(db)
}

#file opening
cell = fread(args$cell, col.names = c("bc","gene_name","transcript_name","umi","frag_prob_weighted"), nThread =1) %>% distinct() %>% na.omit()

#Perform EM algorithm
print("EM algorithm...")
genetabs = split(cell, cell$gene_name)
res = Reduce(rbind, lapply(genetabs, function(x) em_algorithm(x)))

#join with cell information
cell = left_join(cell,res)
cell = cell %>% dplyr::distinct(bc,gene_name,transcript_name,rel_abund)

#rounding to 3dec threshold
cell$rel_abund = round(cell$rel_abund, 3)

#writing
fwrite(cell,file = args$output_path, sep="\t", col.names = F, row.names = F, nThread=1)
