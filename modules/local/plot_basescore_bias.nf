process PLOT_BASESCORE_BIAS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v8':'kubran/odcf_snvcalling:v8' }"

    input:
    tuple val(meta), path(vcf), path(reference_allele_base_qualities), path(alternative_allele_base_qualities)
    val(pdfname)
    val(title)

    output:
    path "*.pdf"           , emit: plot  
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tripletBased_BQRatio_plotter.R \\
        -v $vcf \\
        -r $reference_allele_base_qualities \\
        -a $alternative_allele_base_qualities \\
        -t ${params.basequal} \\
        -o ${prefix}_${pdfname}.pdf \\
        -p Differences \\
        -d "${prefix}${title}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
