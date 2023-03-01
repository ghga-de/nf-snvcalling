// This process is not tested!
process BEDTOOLS_SUBTRACT {
    tag "$meta.id"
    label 'process_medium'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"   
    input:
    tuple val(meta), file(orjinal_somatic_vcf), file(filtered)

    output:
    tuple val(meta), path("*_to_10_removedByMedian*.vcf")        , emit: subtracted_file      
    path  "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def suffix     = params.rerunfiltering ? "1": ""
    def outfile    = "${prefix}_somatic_functional_snvs_conf_${params.min_confidence_score}_to_10_removedByMedian${params.median_filter_threshold}Filter.vcf"
    
    """
    mv ${orjinal_somatic_vcf} ${orjinal_somatic_vcf}.bed
    mv ${filtered} ${filtered}.bed
    touch $outfile
    grep '^#' ${orjinal_somatic_vcf}.bed >$outfile; bedtools subtract -a ${orjinal_somatic_vcf}.bed -b ${filtered}.bed >>$outfile  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}