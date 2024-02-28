

/* Internal priming filtering processing:
    set of process for processing of internal priming filtering ops
*/


process ip_splitting {
    tag "${chr}"
	publishDir "./results/internalp_filtering/ipdb_splitted"
    maxForks params.cpus
    cpus params.cpus
    cache true

    input:
        tuple val(chr), path(ipref)

    output:
	    tuple val(chr), path("${chr}_sp.ipdb")
	
    script:
        """
        #select chromosome
        gawk '{ if(\$1=="${chr}") print }' ${ipref} > ${chr}_sp.ipdb
        """
}


process ip_filtering {
    tag "${sample_id}, ${chr}, ${ipdb}"
	publishDir "./results/internalp_filtering/internalp_filtered/${sample_id}"
    maxForks params.cpus
    cpus params.cpus
    cache true

    input:
        tuple val(sample_id), val(chr), path(reads), path(ipdb)

    output:
	    tuple val(sample_id), val(chr), path("${chr}_mapped.ipdb"), path("${chr}_unique.reads"), path("${chr}.readid"), path("${chr}_ipf.reads")
	
    script:
        """
        #filtering
        Rscript ${baseDir}/src/ip_filtering.R ${reads} ${ipdb} ${params.isoform_end_ip_threshold} ${chr}.readid ${chr}_unique.reads ${chr}_mapped.ipdb ${chr}_ipf.reads
        """
}