process PER_CHROM_PLOT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(vcf)
    val(exy)

    output:
    tuple val(meta), path('*.txt.tmp')   , emit: distance_tmp        
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    mutationDistance.py \\
        --inf=$vcf \\
        --outf=${prefix}_somatic_mutation_dist_conf_${params.min_confidence_score}_to_10.txt.tmp \\
        --alleleFreq=$params.allele_freq \\
        $exy

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2.7 --version | sed 's/Python //g')
    END_VERSIONS
    """
}