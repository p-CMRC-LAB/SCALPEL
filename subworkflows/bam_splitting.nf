
/* - Splitting of the BAM file by chromosome ID and Processing of the BAMs:
===========================================================================
*/

process bam_splitting {
        tag "${sample_id}, ${chr}, ${bc_path}"
        publishDir "./results/reads_processing/bam_splitting/${sample_id}"
        cache true
        label 'small_mem'

	input: 
            tuple val(chr), val(sample_id), path(bam), path(bai), path(dge_matrix), val(bc_path)
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
                    samtools view --subsample ${params.subsample} -b ${bam} ${chr} -D XC:${bc_path} --keep-tag "XC,XM" | samtools sort > tmp.bam
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
                    samtools view --subsample ${params.subsample} -b ${bam} ${chr} --keep-tag "XC,XM" | samtools sort > tmp.bam
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
                    samtools view --subsample ${params.subsample} -b ${bam} ${chr} -D CB:${bc_path} --keep-tag "CB,UB" | samtools sort > tmp.bam
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
                    samtools view --subsample ${params.subsample} -b ${bam} ${chr} --keep-tag "CB,UB" | samtools sort > tmp.bam
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
