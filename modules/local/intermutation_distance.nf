process INTERMUTATION_DISTANCE {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(vcf)
    each path(chr_file)
    val(chr_prefix)

    output:
    tuple val(meta), path('*.pdf')   , emit: plot        
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    intermutationDistance_Coord_color.r \\
        -i $vcf \\
        -s $prefix \\
        -o snvs_${prefix}_intermutation_distance_conf_${params.min_confidence_score}_to_10.pdf \\
        -a "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \\
        -l $chr_file \\
        -p "${chr_prefix}" \\
        -u ""

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}