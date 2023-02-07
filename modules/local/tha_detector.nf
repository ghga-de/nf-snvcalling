process THA_DETECTOR {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(somatic_vcf)

    output:
    tuple val(meta), path('*.txt')        , emit: tha_file        
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    THA_SCORE=`determine_THA_score.R  -i ${somatic_vcf}`
    [[ \$(echo "${THA_SCORE} > ${params.tha_score_threashold}" | bc -l) ]] && echo -e "THA score\\t\${THA_SCORE}\\n" >${prefix}_is_THA_affected.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}