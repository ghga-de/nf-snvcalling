process MAKE_MAF_INPUT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"
    
    input:
    tuple val(meta), file(vcf)

    output:
    tuple val(meta), path('*.txt')   , emit: maf_values        
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    makeMAFinput.pl \\
        $vcf \\
        $params.min_cov \\
        $params.min_confidence_score > ${prefix}_MAF_conf_${params.min_confidence_score}_to_10.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}