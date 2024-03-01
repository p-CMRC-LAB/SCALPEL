


process extraction_readids {
    tag "${sample_id}, ${chr}"
    cache true

    input:
        tuple val(sample_id), val(chr), path(reads)
    output:
        tuple val(sample_id), path("${chr}_reads_selected.txt")

    script:
        """
        awk -v OFS="\\t" 'NR > 1 {print \$14,\$19,\$21}' ${reads} | sort -u -T . > ${chr}_reads_selected.txt
        """
}


process merge_extracted_readsids {
    tag "${sample_id}"
    publishDir "./results/final_results/${sample_id}"
    cache true

    input:
        tuple val(sample_id), path(reads)
    output:
        tuple val(sample_id), path("${sample_id}_reads_selected.RDS")

    script:
        """
        #!/usr/bin/env Rscript
        library(data.table)
        library(dplyr)

        #get files
        all.files = list.files(".", pattern="*_reads_selected.txt", full.names=T)
        #reading
        all.reads = lapply(all.files, function(x) fread(x, col.names = c("read.id","gene_name","transcript_name"), nThread=1)) %>% rbindlist()
        #writing
        saveRDS(all.reads, file='${sample_id}_reads_selected.RDS')
        """
}


process bam_filtering {
	tag "${sample_id}, ${chr}"
    cache true
	
    input:
        tuple val(sample_id), val(chr), path(bam), path(ipdb), path(readid)
    output:
    	tuple val(sample_id), path("${chr}_ft.bam"), path(ipdb)

	script:
	"""
	#Filter bam file based on read id selected
	samtools view -@ 1 -b -N ${readid} ${bam}  > ${chr}_ft.bam
	"""
}


process merging_results {
    tag "${sample_id}"
	publishDir "./results/final_results/${sample_id}", mode: 'copy'
    cache true

	input:
        tuple val(sample_id), file(bamfiles), file(ipdbs)
	output:
        tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), path("${sample_id}_internalp.txt")
    
	script:
        """

        # FILTERED BAM FILES
        # ==================

        #merge all the bams files (1)
        #========================
        samtools merge -@ 1 -f -o final1.bam ${bamfiles}
        samtools sort -@ 1 final1.bam > ${sample_id}.filtered.bam
        samtools index -@ 1 ${sample_id}.filtered.bam
        rm final1.bam


        # ANNOTATIONS
        
        #collect all ips (2)
        echo "seqnames	start_ip	end_ip	strand	exon_id	start	end	gene_name	transcript_name	exon_ind	start_rel	end_rel	rel_start_ip	rel_end_ip	dist_END_ip" > ${sample_id}_internalp.txt
        cat *_mapped.ipdb >> ${sample_id}_internalp.txt
        """
}
