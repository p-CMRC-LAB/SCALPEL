

/* Reads_processing:
    set of process for processing of aligned reads and filtering ops
*/


process bam_splitting {
	tag "${sample_id}, ${chr}"
	publishDir "./results/reads_processing/bam_splitting/${sample_id}"
    cache true

	input: 
        tuple val(chr), val(sample_id), file(sample)

	output:
        tuple val(sample_id), val("${chr}"), path("${chr}.bam")

	script:
        if( params.sequencing == "dropseq" )
            if ( params.barcodes != null )
                """
                    #Dropseq
                    #=======
                    #input files
                    #Filter reads , Remove duuplicates and split by chromosome
                    #a) filter reads...
                    samtools view -@ 1 -b ${sample[0]} ${chr} -D XC:${params.barcodes} --keep-tag "XC,XM" | samtools sort > tmp.bam
                    #Remove all PCR duplicates ...
                    samtools markdup tmp.bam ${chr}.bam -r --barcode-tag XC --barcode-tag XM
                    #delete all empty files ...
                    rm tmp.bam
                    samtools view ${chr}.bam | head -2 > check
                    if [ -s check ]; then
                        echo "ok"
                    else
                        rm -f ${chr}.bam
                    fi
                """
            else
                """
                    #Dropseq
                    #=======
                    #input files
                    #Filter reads , Remove duuplicates and split by chromosome
                    #a) filter reads...
                    zcat -f ${sample[2]} > bc.txt
                    samtools view -@ 1 -b ${sample[0]} ${chr} -D XC:bc.txt --keep-tag "XC,XM" | samtools sort > tmp.bam
                    #Remove all PCR duplicates ...
                    samtools markdup tmp.bam ${chr}.bam -r --barcode-tag XC --barcode-tag XM
                    #delete all empty files ...
                    rm tmp.bam
                    samtools view ${chr}.bam | head -2 > check
                    if [ -s check ]; then
                        echo "ok"
                    else
                        rm -f ${chr}.bam
                    fi
                """
        else if( params.sequencing == "chromium" )
            if ( params.barcodes != null )
                """
                    #Chromium_seq
                    #============
                    #Filter reads , Remove duuplicates and split by chromosome
                    samtools view -@ 1 -b ${sample[0]} ${chr} -D CB:${params.barcodes} --keep-tag "CB,UB" | samtools sort > tmp.bam
                    #Remove all PCR duplicates ...
                    samtools markdup tmp.bam ${chr}.bam -r --barcode-tag CB --barcode-tag UB
                    #delete all empty files ...
                    rm tmp.bam
                    samtools view ${chr}.bam | head -2 > check
                    if [ -s check ]; then
                        echo "ok"
                    else
                        rm -f ${chr}.bam
                    fi
                """
            else
                """
                    #Chromium_seq
                    #============
                    #Filter reads , Remove duuplicates and split by chromosome
                    zcat -f ${sample[2]} | sort -u > bc.txt
                    samtools view -@ 1 -b ${sample[0]} ${chr} -D CB:bc.txt --keep-tag "CB,UB" | samtools sort > tmp.bam
                    #Remove all PCR duplicates ...
                    samtools markdup tmp.bam ${chr}.bam -r --barcode-tag CB --barcode-tag UB
                    #delete all empty files ...
                    rm tmp.bam
                    samtools view ${chr}.bam | head -2 > check
                    if [ -s check ]; then
                        echo "ok"
                    else
                        rm -f ${chr}.bam
                    fi
                """
}


process bedfile_conversion{
	tag "${sample_id}, ${chr}"
	publishDir "./results/reads_processing/bedfile_conversion/${sample_id}"
    cache true

	input:
        tuple val(sample_id), val(chr), path(bam)

	output:
        tuple val(sample_id), val(chr), path("${chr}.bed")

	script:
        """
            #Convertion of bam files to bed files (1)
            bam2bed --all-reads --split --do-not-sort < ${bam} | gawk -v OFS="\\t" '{print \$1,\$2,\$3,\$6,\$4,\$14"::"\$15}' > temp_file

            #Extract relevant metadata infos and sort for overlapping purposes (2)
            Rscript ${baseDir}/src/subset_metadata.R temp_file ${chr}.bed

            #delete tmp files (3)
            rm temp_file
        """
}


process reads_mapping_filtering {
	tag "${sample_id}, ${chr}, ${bed}, ${exons}"
	publishDir "./results/reads_processing/mapping_filtering/${sample_id}"
    cache true

	input:
        tuple val(sample_id), val(chr), path(bed), path(exons), path(exons_unique)
	
    output:
        tuple val(sample_id), val(chr), path("${chr}.reads"), path("encoding.txt")

	script:
        """
            Rscript ${baseDir}/src/mapping_filtering.R ${bed} ${exons} ${params.dt_threshold} ${chr}.reads
        """
}






