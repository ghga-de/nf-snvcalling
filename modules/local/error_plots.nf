process ERROR_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"

    input:
    tuple val(meta), path(vcf)
    val(errortype)
    val(pdfname)
    val(errorfilename)
    val(plottitle)

    output:
    tuple val(meta), path("*txt"), emit: error_matrix
    path  "*.pdf"                , emit: plot
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    createErrorPlots.py --vcfFile=$vcf \\
        --referenceFile=NA \\
        --outputFile=${prefix}_${pdfname}.pdf \\
        --errorType=$errortype \\
        --errorFile=${prefix}_${errorfilename}.txt \\
        --plot_title="${plottitle}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
