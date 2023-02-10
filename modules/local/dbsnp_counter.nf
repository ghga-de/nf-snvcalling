process DBSNP_COUNTER {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(vcf)

    output:
    tuple val(meta), path('*.txt')        , emit: indbsnp        
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    in_dbSNPcounter.pl \\
        $vcf \\
        $params.min_confidence_score > ${prefix}_somatic_in_dbSNP_conf_${params.min_confidence_score}_to_10.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}