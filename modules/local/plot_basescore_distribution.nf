process PLOT_BASESCORE_DISTRIBUTION {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v8':'kubran/odcf_snvcalling:v8' }"

    input:
    tuple val(meta), path(vcf), path(reference_allele_base_qualities), path(alternative_allele_base_qualities)
    val(pdfname)
    val(title)

    output:
    path "*pdf"           , emit: plot  
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotBaseScoreDistribution.R \\
        -v $vcf \\
        -r $reference_allele_base_qualities \\
        -a $alternative_allele_base_qualities \\
        -t ${params.basequal} \\
        -o ${prefix}_${pdfname}.pdf \\
        -d "${prefix} ${title}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
