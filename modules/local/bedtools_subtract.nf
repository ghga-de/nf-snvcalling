// This process is not tested!
process BEDTOOLS_SUBTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2':
        'quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"   
    input:
    tuple val(meta), file(somatic_vcf), file(filtered)

    output:
    tuple val(meta), path("*_to_10_removedByMedian*.vcf")        , emit: somatic_functional      
    path  "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def suffix     = params.rerunfiltering ? "1": ""
    
    """
    grep '^#' ${somatic_vcf} >${filtered}; bedtools subtract -a ${somatic_vcf} -b ${filtered} >>${prefix}_somatic_functional_snvs_conf_${params.min_confidence_score}_to_10_removedByMedian${params.median_filter_threshold}Filter.vcf  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\1/')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}