process ERROR_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*_sequencing_error_matrix.txt"), path("*_sequence_error_matrix.txt")  , emit: error_matrix
    path  "versions.yml"                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    createErrorPlots.py --vcfFile=$vcf \\
        --referenceFile=NA \\
        --outputFile=${prefix}_sequencing_specific_error_plot_before_filter.pdf \\
        --errorType=sequencing_specific \\
        --errorFile=${prefix}_sequencing_error_matrix.txt \\
        --plot_title='Sequencing strand bias before guanine oxidation filter'

    createErrorPlots.py --vcfFile=$vcf \\
        --referenceFile=NA \\
        --outputFile=${prefix}_sequence_specific_error_plot_before_filter.pdf \\
        --errorType=sequence_specific \\
        --errorFile=${prefix}_sequence_error_matrix.txt \\
        --plot_title='PCR strand bias before guanine oxidation filter'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
