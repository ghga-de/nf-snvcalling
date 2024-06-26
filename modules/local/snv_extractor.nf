process SNV_EXTRACTOR {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta), path('*_somatic_functional_snvs_conf_*_to_10*')        , emit: somatic_functional
    tuple val(meta), path('*_somatic_snvs_conf_*_to_10*')                   , emit: somatic_snv
    tuple val(meta), path('*_somatic_functional_ncRNA_snvs_conf_*_to_10*')  , emit: somatic_ncrna              
    tuple val(meta), path('*_germline_functional_snvs_conf_*_to_10*')       , emit: germline_functional  
    path  "versions.yml"                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    snv_extractor_v1.pl \\
        --infile=$vcf \\
        --minconf=$params.min_confidence_score \\
        --pid=snvs_$prefix \\
        $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}