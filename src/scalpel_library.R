





# [Script for Function needed in Scalpel analysis]
# ------------------------------------------------
library(GenomicRanges)
library(Gviz)

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
  all_tests = pbapply::pblapply(genes_names, function(x){
    #a- get count table
    a = ALL_expression_by_GENE[[x]][,conds]
    b = apply(a,2, function(x) x/sum(x))
    c = a[names(which(rowSums(b > threshold_tr_abundance) >= 1)),]
    if(nrow(c) > 1){
      d = suppressWarnings(chisq.test(c, correct = T, simulate.p.value = F))
      c$gene = x
      c$p_value = d$p.value
      c$p_value.adjusted = p.adjust(c$p_value,method = 'BH')
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
  
  #filter
  RES_TAB_SIGNIF = RES_TAB %>% filter(p_value.adjusted < pval_adjusted)
  RES_TAB_SIGNIF$gene_tr = rownames(RES_TAB_SIGNIF)
  RES_TAB_SIGNIF$transcript = stringr::str_split_fixed(RES_TAB_SIGNIF$gene_tr,pattern = "\\*\\*\\*",n=2)[,2]
  RES_TAB_SIGNIF = RES_TAB_SIGNIF %>% arrange(p_value.adjusted,gene)
  
  #5/ Return
  return(RES_TAB_SIGNIF)
}







CoveragePlot = function(genome_gr=NULL, gene=NULL, bamfiles=NULL, bamnames=NULL,
                        seuratobj=NULL, condition=NULL, iptab = NULL, patab = NULL,
                        genome=NULL, transcriptAnnotation="name", transcripts_tofilter=NULL,
                        samtoolsbin = "samtools") {
    
    # Function to plot Transcript Coverage
    # #***********************************
    # 
    # required args:
    #     - genome_gr [Granges object]
    #     - gene [str]
    #     - bamfiles [str]
    #     - bamnames [str]
    #     - seuratobj
    # 
    # optional args:
    #     -
    #     -
    #     -
    
    
    # Check input args
    # ----------------
    print('cheking inputs args...')
    if(is.null(genome_gr) || is.null(gene) || is.null(bamfiles) || is.null(bamnames) || is.null(seuratobj)){
        stop("Error in input args provided, check required args (genome_gr, gene, bamfiles, bamnames, seuratobj) !!")
    }
    
    if(!is.null(condition)){
        if(length(table(seuratobj[[condition]])) != length(bamfiles)){
            stop("Error: The number of provided bamfiles does not match the level of condition specified in the seurat object")
        }else{
            size_samples = table(seuratobj[[condition]])[bamnames]
            print("sample sizes...")
            print(size_samples)
        }
    }else{
        size_samples = rep(1, length(bamfiles))
    }
    
    if(is.null(genome)){
        stop("Please: Specify genome (hg39 or mm10)")
    }
    

    # Filter Grange genome table and get gene metada (1)
    # ----------------------------------------------
    print('Filter Grange genome table and get gene metada...')
    genein = genome_gr[(genome_gr$gene_name==gene) & (genome_gr$type=="exon")]
    if(length(genein)==0){stop("Error: Empty gene table.... Gene provided not found into the GenomeGrange object! Check gene provided or GTF!")}
    if(!(is.null(transcripts_tofilter))){genein = genein[genein$transcript_name %in% transcripts_tofilter]}
    transcripts = GenomicFeatures::transcripts(GenomicFeatures::makeTxDbFromGRanges(genein))
    
    #define sets
    chrom = as.character(seqnames(genein)@values)
    start = min(start(genein))
    end = max(end(genein))
    gene_strand = as.character(strand(genein)@values)
    
    #Internal priming track
    if(is.null(iptab)){
        atrack = AnnotationTrack(name = "ip_sites", background.title="cornflowerblue", col="red", background.title="cornflowerblue", chromosome = chrom, genome = genome)
    }else{
        print("internal priming table provided...")
        ip_gene = iptab %>% filter(gene_name == gene) %>% dplyr::select(seqnames, start_ip, end_ip, strand) %>% distinct()
        ip_gene = ip_gene %>% filter(start_ip >= start & end_ip <= end)
        colnames(ip_gene) = c("seqname","start","end","strand")
        ip_gene$start = as.numeric(ip_gene$start)
        ip_gene$end = as.numeric(ip_gene$end)
        if(nrow(ip_gene)>0){
            print(paste0(gene, ": associated ip_pos founded..."))
            atrack = AnnotationTrack(ip_gene, name = "ip_sites", chromosome = chrom, col="red", background.title="cornflowerblue")
        }else{
            print(paste0(gene, ": no ip pos"))
            atrack = AnnotationTrack( name = "ip_sites", col="red", background.title="cornflowerblue")   
        }
    }
    
    
    #PolyA sites track
    if(is.null(patab)){
        patrack = AnnotationTrack(name = "polyA_sites", background.title="cornflowerblue", chromosome = chrom, genome = genome)
    }else{
        print("PA sites table provided...")
        start_lim = start
        end_lim = end
        patab = patab %>% dplyr::filter((type=="polyA_site") & (start >= start_lim & end <= end_lim)) %>%
            dplyr::select(seqnames, start, end, strand, type) %>% distinct()
        patrack = AnnotationTrack(patab, name = "polyA_sites", chromosome = chrom, col="blueviolet", background.title="cornflowerblue")
    }

    
    # Create Region Track object (2)
    # --------------------------
    print("Create Region Track object...")
    genein.tb = data.frame(genein)
    
    #check columns names issues
    colnames(genein.tb) = str_replace(colnames(genein.tb),"gene_biotype","gene_type")
    colnames(genein.tb) = str_replace(colnames(genein.tb),"transcript_biotype","transcript_type")
    
    
    #selection
    genein.tb = genein.tb[,c("seqnames","start","end","width","strand","gene_type","gene_id",
    "exon_id","transcript_id","gene_name","transcript_name")]
    fwrite(genein.tb[,c("seqnames","start","end","strand")] %>% arrange(start,end), file = "exon.bed", sep="\t")
    #rename
    colnames(genein.tb) = c("Chromosome","start","end","width","strand","feature","gene","exon",
    "transcript","gene_name","transcript_name")
    if(transcriptAnnotation == "name"){
    print("transcriptAnnotation [transcriptName]")
    genein.tb$transcript = genein.tb$transcript_name
    }else{
    print("transcriptAnnotation [transcriptID]")
    }
    #Track
    gtrack = GeneRegionTrack(genein.tb, chromosome = chrom, genome = "hg38",
    name = gene, transcriptAnnotation = "transcript",
    background.title = "darkblue", just.group="below")

    
    # Create Coverate Track objects (3)
    # -----------------------------
    print("Create Coverate Track objects...")
    tryCatch(
        expr = {
            DTRACKS = lapply(1:length(bamfiles), function(x){
                print(bamfiles[x])
                #open bamfile
                bamfile = Rsamtools::BamFile(bamfiles[x])
                #get coverage
                coveragein = GenomicRanges::coverage(bamfile, param=Rsamtools::ScanBamParam(which = transcripts))
                #define Granges
                gr = GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges::IRanges(start:end, width = 1), strand = gene_strand)
                mcols(gr) = (as.numeric(coveragein[[chrom]])[start:end])# / size_samples[x])
                #return Dtrack object
                RET = list()
                RET$OBJ = Gviz::DataTrack(gr, type = "l", genome = genome, name = bamnames[x], fill.histogram="black", col.histogram = "darkgrey", 
                                          background.title = "black")
                RET$YMAX = max(mcols(gr)$X)
                
                return(RET)
            })
            YMAX = max(unlist(lapply(DTRACKS, function(x) x$YMAX)))
            DTRACKS = unlist(lapply(DTRACKS, function(x) x$OBJ))
            
            # [Chromosome Track]z
            ctrack = GenomeAxisTrack()

            # [Final Plot]
            # ------------
            print("Plotting...")
            plotTracks(c(ctrack,atrack, patrack, gtrack, DTRACKS), ylim = c(0,YMAX))
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
    })
}



plot_relativeExp = function(seurat.obj, features, group.var, levels.group=NULL, ...){
    #Function to plot with geom_boxplot relative isoform expression
    Reduce(`+`, lapply(features, function(x) {
        subsetted_expression = seurat.obj[x,]@meta.data %>% select(nCount_RNA, `group.var`)
        subsetted_expression$nCount_RNA = (subsetted_expression$nCount_RNA / sum(subsetted_expression$nCount_RNA)) * 100
        if(is.null(levels.group)){
            print("no level reordering")
        }else{
            subsetted_expression[,group.var] = factor(subsetted_expression[,group.var], levels = levels.group)
        }
        ggplot(subsetted_expression, aes_string(group.var, "nCount_RNA")) +
            geom_boxplot(aes_string(fill=group.var)) +
            ggtitle(label = x) + theme_classic(base_size = 13)
    })) + plot_layout(ncol = length(features))
}



CoverPlot = function(genome_gr, gene_in, genome_sp, bamfiles, distZOOM=NULL, annot_tab=NULL,
                     transcripts_in=NULL, filter_trs=F, samtools.bin="samtools"){
  #check args
  print(genome_sp)
  if(is.null(genome_sp) || is.null(genome_gr) || is.null(bamfiles)){
    stop("Error ! Check input args : genome.sp / genome.gr / bamfiles")
  }

  #0. Ideogram
  axisTrack <- GenomeAxisTrack(genome=genome_sp)

  #1. Build genome table & Track
  print("GenomeTrack building...")
  gtab.gr = genome_gr[genome_gr$gene_name==gene_in]
  gtab.gr$colors = "orange"
  #filter transcripts provided in case
  if(filter_trs){
    print("filtering...")
    print(transcripts_in)
    if(is.null(transcripts_in)){stop("Error! Provide list of transcripts: transcripts_in")}
    gtab.gr = gtab.gr[gtab.gr$transcript_name %in% transcripts_in]
  }
  gtab.gr$colors = as.character(gtab.gr$colors)
  gtab.gr$transcript_id = paste(gtab.gr$transcript_id, gtab.gr$transcript_name,sep=" - ")
  gtab = data.frame(gtab.gr) %>%
    dplyr::filter(type %in% c("UTR","CDS","exon")) %>%
    dplyr::distinct(seqnames,start,end,width,strand,gene_type,gene_id,
                    exon_id,transcript_id,gene_name,transcript_name,type) %>%
    dplyr::rename(Chromosome="seqnames",feature="gene_type",
                  gene="gene_id",exon="exon_id",transcript="transcript_id") %>%
    group_by(transcript) %>%
    mutate(check = ifelse("UTR" %in% type, T,F)) %>%
    filter(!(type == "exon" & check == T)) %>%
    ungroup()
  gtab$feature = gtab$type
  gtab$feature = stringr::str_replace(gtab$feature,"UTR","utr")
  gtab$feature = stringr::str_replace(gtab$feature, "exon", "utr")
  gtab$colors = "orange"
  if(!is.null(transcripts_in)){
    print("highliting transcripts...")
    print(transcripts_in)
    # gtab$colors[which(gtab$transcript_name %in% transcripts_in)] = "blue"
    gtab = arrange(gtab, transcript, start, desc(end))
    # gtab = Reduce(rbind, lapply(split(gtab, gtab$transcript), function(x) x[order(x$start),]))
  }
  #define sets
  chrom = gtab$Chromosome[1]
  starts = min(gtab$start)
  ends = max(gtab$end)
  strands = as.character(gtab$strand[1])

  #track
  print(gtab)
  dplyr::select(gtab, Chromosome,start,end) %>%
    mutate(start = start - 25, end = end + 25) %>%
  data.table::fwrite(file="./coords.txt", sep="\t", row.names = F, col.names = F)
  Gtrack = Gviz::GeneRegionTrack(GenomicRanges::makeGRangesFromDataFrame(gtab, keep.extra.columns = T), chromosome = chrom,
                                 name = gene_in, transcriptAnnotation = "transcript", 
                                 just.group="below", genome=genome_sp, fill="orange", color="black", col = "black", 
                                 background.title="darkmagenta", fontsize.group=16, 
                                 fill=gtab$colors)

  #1bis : Annotation Track building
  if(!is.null(annot_tab)){
    curr.tab = dplyr::filter(annot_tab, start>=starts & end<=ends)
    Atrack = AnnotationTrack(GenomicRanges::makeGRangesFromDataFrame(curr.tab), name = "Peaks",
                             chromosome = chrom, fill="olivedrab", id=curr.tab$name,
                             background.title="black", featureAnnotation = "id", fontcolor.group="black")
  }else{
    Atrack = AnnotationTrack(chromosome = chrom, background.title="black", name = "Peaks")
  }

  #2. Build Alignment Track
  print("DataTack building...")
  print("oooook")
  Atracks.res = lapply(1:length(bamfiles), function(x){
    print(names(bamfiles)[x])
    #get coverage
    system(paste0(samtools.bin, " view -b --region-file coords.txt ", bamfiles[x], " > current.bam"))
    cov.exp = system(paste0(samtools.bin," depth -b coords.txt current.bam > current.cov"))
    cov.tab = fread("current.cov", col.names = c("seqnames","start","depth")) %>% dplyr::filter(depth>=0)
    #dataTrack
    curr.track = DataTrack(start = cov.tab$start, width=1, data = cov.tab$depth, chromosome = chrom, genome=genome_sp,
                           type=c("hist"), background.title="coral4", name = names(bamfiles)[x], col.histogram="blue")
    return(list(max.depth=max(cov.tab$depth), track=curr.track))
  })
  Atracks = lapply(Atracks.res, function(x) x[[2]])
  YMAX = max(unlist(lapply(Atracks.res, function(x) x[[1]])))
  system("rm current.cov")
  system("rm current.bam")
  print(Atracks)

  #plotting...
  print("Plotting...")
  in.tracks = c(axisTrack, Atrack,Gtrack, Atracks)

  #in case of ZOOMING option
  size.tracks = c(0.3,0.1,0.6,rep(0.9/length(bamfiles), length(bamfiles)))
  if(is.null(distZOOM)){
    plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX))
  }else{
    if(strands == "+"){
      size.tracks = c(0.1,0.4,0.6,rep(0.5/length(bamfiles), length(bamfiles)))
      print(size.tracks)
      plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX), from=ends-distZOOM, just.group="right")}
    if(strands == "-"){
      size.tracks = c(0.1,0.4,0.6,rep(0.5/length(bamfiles), length(bamfiles)))
      print(size.tracks)
      plotTracks(in.tracks, sizes = size.tracks, ylim = c(0,YMAX), to=starts+distZOOM, just.group="left")}
  }
}




