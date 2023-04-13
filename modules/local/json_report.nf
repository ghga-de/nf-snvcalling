process JSON_REPORT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"
    
    input:
    tuple val(meta), file(somatic_vcf), file(insnp_file)

    output:
    tuple val(meta), path('*is_THA_affected.txt')   , emit: tha
    tuple val(meta), path('*.json')                 , emit: json
    tuple val(meta), path('*.pdf')                  , emit: plot        
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"

    """
    final_plots_and_json.sh \\
        -p $prefix \\
        -i $somatic_vcf \\
        -s $insnp_file \\
        -t $params.min_cov \\
        -v $params.min_confidence_score \\
        -b $params.tha_score_threshold \\
        -r ''

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}