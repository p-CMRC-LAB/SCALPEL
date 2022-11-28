

#Script for visualization of genome coverage
#Author: Franz AKE
#July 2022
genome_cover = function(genome_gr, bamfiles, bamnames, gene_in, targets_trs=c(), sample_sizes=c(),
                        extra_left=0, extra_right=0){
  #Paths
  my_Rlibs = "~/CEPH/R_PACKAGES/"
  
  #Create a Gene track (2)
  #*******************
  gene_gr = genome_gr[genome_gr$gene_name == gene_in]
  #Merge transcript_id and transcript_name
  gene_gr$transcript_id = paste(gene_gr$transcript_name, gene_gr$transcript_id, sep="/")
  #filter specific transcripts
  if(length(targets_trs)!=0){
    gene_gr = gene_gr[gene_gr$transcript_name %in% targets_trs]
  }
  #get coordinates
  start_gene = min(start(gene_gr))
  end_gene = max(end(gene_gr))
  strand_gene = unique(as.character(strand(gene_gr)))
  seqlevels(gene_gr) = na.omit(as.character(seqnames(gene_gr)))[1]
  #Make TxDb from Granges objects 
  gene_txdb = suppressWarnings(makeTxDbFromGRanges(gene_gr))
  if(strand_gene == "-"){
    gtrack = GeneRegionTrack(gene_txdb, transcriptAnnotation = 'symbol', name = gene_in, just.group = 'above', fontcolor.group="blue",
                             arrowHeadMaxWidth= 5, cex.group=1, fontsize = 6, fontsize.group = 6, background.title = "cornflowerblue", size=c(30))
  }else{
    gtrack = GeneRegionTrack(gene_txdb, transcriptAnnotation = 'symbol', name = gene_in, just.group = 'above', fontcolor.group="blue",
                             arrowHeadMaxWidth= 5, cex.group=1, fontsize = 6, fontsize.group = 6, background.title = "cornflowerblue", size = c(30))
  }
  ideoTrack <- IdeogramTrack(genome = "hg38", chromosome = seqlevels(gene_gr), from=start_gene, to=end_gene, cex=1)
  
  # Get Genome read coverage (3)
  #*************************
  param = Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(seqnames = seqlevels(gene_gr),
                                                                 ranges = IRanges::IRanges(start_gene - 0, width = end_gene - start_gene + 0), strand = strand_gene))
  dtracks_res = lapply(1:length(bamfiles), function(x){
    print(bamfiles[x])
    cov = GenomicAlignments::coverage(bamfiles[x], param = param)[[seqlevels(gene_gr)]]
    dtrack = DataTrack(start=start(cov), end = end(cov), data = (runValue(cov)/sample_sizes[x]) * 1, genome="hg38", type=c('h'), chromosome = seqlevels(gene_gr),
                       cex = 3, cex.axis=0.5,cex.main=1, cex.title = 1.3, name = bamnames[x], background.title = "gray", size=c(80))
    max_val = max(runValue(cov)/sample_sizes[x]) * 1
    return(list(max_val, dtrack))
  })
  dtracks = unlist(lapply(dtracks_res, function(x) x[2]))
  max_y = max(unlist(lapply(dtracks_res, function(x) x[1])))
  
  #Plotting (4)
  #********
  options(repr.plot.width=40, repr.plot.height=20)
  axisTrack = GenomeAxisTrack(range=IRanges(start = c(start_gene), end = c(end_gene), names = c("Gene_range")), cex=1, exponent=4)
  if(strand_gene == "+"){
    plotTracks(c(axisTrack, dtracks, gtrack), from=start_gene-extra_left, to=end_gene+2e3+extra_right,
               background.panel = "#FFFEDB", background.title = "darkblue", showId=T, add53 = TRUE, add35 = TRUE, littleTicks = TRUE)    
  }else{
    plotTracks(c(axisTrack, dtracks, gtrack), from=start_gene-2e3-extra_left,
               to=end_gene+extra_right, showId=T, add53 = TRUE, add35 = TRUE, littleTicks = TRUE,
               background.panel = "#FFFEDB", background.title = "darkblue")
  }
}
  
  
  
  
  
  
  
