process MUTATION_DISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.txt')   , emit: distance      
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    mutationDistance.py \\
        --inf=$vcf \\
        --outf=snvs_${prefix}_somatic_mutation_dist_conf_${params.min_confidence_score}_to_10.txt \\
        --alleleFreq=$params.allele_freq 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}