process JSON_REPORT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(vcf), file(index)

    output:
    tuple val(meta), path('snv_*_somatic_functional_snvs_conf_*_to_10.vcf')        , emit: somatic_functional
    tuple val(meta), path('snv_*_somatic_snvs_conf_*_to_10.vcf')                   , emit: somatic_snv
    tuple val(meta), path('snv_*_somatic_functional_ncRNA_snvs_conf_*_to_10.vcf')               
    tuple val(meta), path('snv_*_germline_functional_snvs_conf_*_to_10.vcf')         
    path  "versions.yml"                                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    snv_extractor_v1.pl \\
        --infile=$vcf \\
        --minconf=$params.min_confidence_score \\
        --pid=snv_$meta.id \\
        $args 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}