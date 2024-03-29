process ERROR_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(vcf)
    val(errortype)
    val(pdfname)
    val(errorfilename)
    val(plottitle)

    output:
    tuple val(meta), path("*.txt") ,  emit: error_matrix
    tuple val(meta), path("*.pdf") ,  emit: plot
    path  "versions.yml"           ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    createErrorPlots.py --vcfFile=$vcf \\
        --referenceFile=NA \\
        --outputFile=snvs_${prefix}_${pdfname}.pdf \\
        --errorType=$errortype \\
        --errorFile=snvs_${prefix}_${errorfilename}.txt \\
        --plot_title="${plottitle}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
