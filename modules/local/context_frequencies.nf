process CONTEXT_FREQUENCIES {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"
    
    input:
    tuple val(meta), file(vcf)

    output:
    tuple val(meta), path('*.txt')        , emit: seq_contex        
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    SNV_context_frequencies.pl \\
        $vcf \\
        $params.min_confidence_score > ${prefix}_snvs_with_context_conf_${params.min_confidence_score}_to_10.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\1/')
    END_VERSIONS
    """
}