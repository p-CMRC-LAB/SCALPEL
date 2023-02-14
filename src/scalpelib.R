



# [Script for Function needed in Scalpel analysis]
# ------------------------------------------------


Find_isoforms = function(seurat.obj, pval_adjusted=0.05, condition="orig.ident", assay="RNA", threshold_tr_abundance = 0.15){
  # ---------------------------------------------------------------------------
  #Function to Find differentially expressed isoforms in the conditions defined
  # ---------------------------------------------------------------------------

  print("processing...")

  #1/ Find genes with at least 2 transcripts
  print("Find genes with at least 2 transcripts...")
  genes_tr_tab = (rownames(seurat.obj) %>% stringr::str_split_fixed(pattern = "\\*\\*\\*", n = 2)) %>% data.table()
  genes_tr_tab$gene_tr = rownames(seurat.obj)
  colnames(genes_tr_tab) = c("gene", "transcript", "gene_transcript")
  counts_genes_tab = genes_tr_tab$gene %>% table() %>% data.table() %>% dplyr::filter(N > 1)
  genes_tr_tab_filtered = genes_tr_tab %>% dplyr::filter(gene %in% counts_genes_tab$.) %>% arrange(gene_transcript)
  
  #2/ Get Isoforms expression in the condition defined
  print("Get isoforms expression in the condition defined...")
  ALL_expression = AggregateExpression(seurat.obj, features = genes_tr_tab_filtered$gene_transcript,
                                       assays = assay, group.by = condition, verbose = T,
                                       slot = 'count')[[assay]] %>% data.frame()
  ALL_expression$only_gene = (rownames(ALL_expression) %>% stringr::str_split_fixed(pattern = "\\*\\*\\*", n = 2))[,1]
  ALL_expression$gene_tr = rownames(ALL_expression)
  #Split all the table by genes
  ALL_expression_by_GENE = split(ALL_expression, ALL_expression$only_gene)
  
  #3/ Perform Chi2 test
  conds = colnames(ALL_expression_by_GENE[[1]])
  conds = conds[-length(conds):(-length(conds)+1)]
  print("Perform Chi2 test on selected genes...")
  genes_names = names(ALL_expression_by_GENE)
  all_tests = lapply(genes_names, function(x){
    #a- get count table
    a = ALL_expression_by_GENE[[x]][,conds]
    b = apply(a,2, function(x) x/sum(x))
    c = a[names(which(rowSums(b > threshold_tr_abundance) >= 1)),]
    if(nrow(c) > 1){
      d = chisq.test(c)
      c$gene = x
      c$p_value = d$p.value
      return(list(c, d))
    }else{
      return(NULL)
    }
  })
  #deleting NULL occurences
  all_tests = all_tests[!sapply(all_tests,is.null)]
  
  #4/ Get tables extracted
  all_tests = lapply(all_tests, function(x) x[[1]])
  RES_TAB = do.call(rbind, all_tests)
  #adjust _pvalue
  print("P.value adjusting......")
  RES_TAB$p_value.adjusted = p.adjust(RES_TAB$p_value,method = 'BH')
  #filter
  RES_TAB_SIGNIF = RES_TAB %>% filter(p_value.adjusted < pval_adjusted)
  RES_TAB_SIGNIF$gene_tr = rownames(RES_TAB_SIGNIF)
  RES_TAB_SIGNIF$transcript = stringr::str_split_fixed(RES_TAB_SIGNIF$gene_tr,pattern = "\\*\\*\\*",n=2)[,2]
  RES_TAB_SIGNIF = RES_TAB_SIGNIF %>% arrange(p_value.adjusted,gene)
  
  #5/ Return
  return(RES_TAB_SIGNIF)
}




CoveragePlot = function(genome, gene_in=NULL, transcripts_in=NULL, bamfiles=NULL, bamnames=NULL,
                        transcritps_showed = "scalpel", seurat.obj=NULL, genome_gr="hg38",
                        condition="orig.ident", left_space=500, rigth_space=500, tr_position="below",
                        cex.tr_name=1.5, scale = TRUE, type_data=c("coverage"), check_ip_pos = FALSE, scalpel_results=NULL){
  # --------------------------------
  # Function to plot isoforms coverage
  # --------------------------------
  
  #check presence of required inputs
  if(is.null(gene_in)){
    return("Enter gene [gene_in]")
  }
  if(is.null(bamfiles) || is.null(bamnames)){
    return("Enter bamfiles path and the names associated [bamfiles / bamnames]")
  }

  if(check_ip_pos==TRUE & is.null(scalpel_results)){
    return("For checking ip position, scalpel results folder path is required")
  }  
  
  #1/ Preprocess genome tab
  #filter exons
  genome = genome[genome$type == "exon",]
  genome = genome[,c("seqnames","start","end","width","strand","gene_type","gene_id",
                     "exon_id","transcript_id","gene_name","transcript_name")]
  #rename
  colnames(genome) = c("Chromosome","start","end","width","strand","feature","gene",
                       "exon","transcript","gene_name","transcript_name")
  #extract gene_tab
  print("extract gene...")
  gene_tab = genome[genome$gene_name==gene_in,]
  #filter or not the transcript processed in scalpel
  if(transcritps_showed=="scalpel"){
    obj.tab = data.table(stringr::str_split_fixed(rownames(seurat.obj), pattern = "\\*\\*\\*", n=2))
    colnames(obj.tab) = c("gene","tr")
    gene_tab = gene_tab %>% filter(transcript_name %in% obj.tab$tr)
  }else if(transcritps_showed=="all"){
    print("ok")
  }
  #Filter transcript specifics
  if(is.null(transcripts_in)){
    print("No specific filtering of input transcripts...")
  }else{
    gene_tab = gene_tab %>% filter(transcript_name %in% transcripts_in)
  }
  gene_tab$symbol = gene_tab$transcript_name
  gene_tab =  gene_tab %>% mutate_if(is.factor,as.character)
  CHR = gene_tab$Chromosome[1]
  start_gene = min(gene_tab$start)
  end_gene = max(gene_tab$end)
  strand_gene = gene_tab$strand[1]

  
  #2/ Create GeneTrack for Gviz
  print("create GeneTrack viz...")
  gtrack = GeneRegionTrack(gene_tab, chromosome = gene_tab$Chromosome[1], genome = genome_gr,
                           name = as.character(gene_tab$gene_name[1]), transcriptAnnotation = "symbol",
                           background.title = "darkblue", cex.group=cex.tr_name, just.group=tr_position, windowSize=5, shape = "arrow")


  #Bis Create Annotation track if ip db
  atrack = AnnotationTrack(name = "IP", background.title = "cornflowerblue")
  if(check_ip_pos==TRUE & (!(is.null(scalpel_results)))){
    #read ip table
    iptab = fread(paste0(scalpel_results,"/reads/ip_filtered/",CHR,".ipp"),sep="\t")
    #filter based on coordinates and strand
    iptab = iptab %>% dplyr::filter(Start_ip>min(gene_tab$start) & Start_ip<max(gene_tab$end))
    if(nrow(iptab)==0){
      atrack = AnnotationTrack(, name = "IP", background.title = "cornflowerblue")
    }else{
      atrack = AnnotationTrack(start = iptab$Start_ip, end = iptab$End_ip, chromosome = CHR, strand = rep(strand_gene,nrow(iptab)), genone = genome_gr, name = "IP", background.title = "cornflowerblue")
    }
  }

  
  #3/ Build coverage
  tryCatch(
    expr = {
      print("create DataTrack viz...")
      CHR = gene_tab$Chromosome[1]
      sample_sizes = table(seurat.obj[[condition]])
      start_gene = min(gene_tab$start)
      end_gene = max(gene_tab$end)
      strand_gene = gene_tab$strand[1]
      param = Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(seqnames = CHR, ranges = IRanges::IRanges(start_gene, width = end_gene - start_gene + 0), strand = strand_gene))
      dtracks_res = lapply(1:length(bamfiles), function(x){
        print(bamfiles[x])
        #readBAM file
        covrg = GenomicAlignments::coverage(bamfiles[x], param=param)[[CHR]]
        coverage_reads = (runValue(covrg))
        # dtrack = DataTrack(start=start(covrg), end = end(covrg), data = coverage_reads, type=c('h'), genome = genome_gr, chromosome = CHR,
        # 					name = bamnames[x], cex.legend=10, cex=10, fontsize.legend=10)
        max_val = max(coverage_reads)
      	dtrack_in = AlignmentsTrack(bamfiles[x], chromosome=CHR, genome=genome_gr, name=bamnames[x], stacking="pack", background.title = "black")
        return(list(max_val, dtrack_in))
      })
      dtracks = unlist(lapply(dtracks_res, function(x) x[2]))
      max_y = max(unlist(lapply(dtracks_res, function(x) x[1])))
      
      #Plotting
      options(ucscChromosomeNames=FALSE)
      #plots
      EXTRA.L = left_space
      EXTRA.R = rigth_space
      
      #tracks
      if (scale == TRUE){
        plotTracks(c(dtracks, atrack, gtrack), from = min(gene_tab$start) - EXTRA.L, to = max(gene_tab$end) + EXTRA.R, type = type_data, ylim=c(0,max_y), 
          sizes = c(rep(0.3, length(dtracks)),0.1,0.5))
        }else{
          plotTracks(c(dtracks, atrack, gtrack), from = min(gene_tab$start) - EXTRA.L, to = max(gene_tab$end) + EXTRA.R, type = type_data, 
          sizes = c(rep(0.5, length(dtracks)),0.1,2))
        }
    },
    error = function(e){
      message("Error ! Use correct inputs !")
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    },
    finally = {
      message('All done, quitting.')
    }
  )
}
