// THA calculation depends in on PV4, DP4, RBS and VDB fileds in the vcf file. bcftool version may change those!
// THA calculation may end with NA
process JSON_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(somatic_vcf), path(insnp_file)

    output:
    tuple val(meta), path('*.txt')                  , emit: txt
    tuple val(meta), path('*.json')                 , emit: json, optional: true
    tuple val(meta), path('*.pdf')                  , emit: plot        
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"

    """
    final_plots_and_json.sh \\
        -p $prefix \\
        -i $somatic_vcf \\
        -s $insnp_file \\
        -t $params.min_cov \\
        -v $params.min_confidence_score \\
        -b $params.tha_score_threshold \\
        -r '' \\
        -j $params.report_json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}