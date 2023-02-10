process MERGE_PLOTS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v9':'kubran/odcf_snvcalling:v9' }"

    input:
    tuple val(meta), file(plot1), file(plot2), file(plot3), file(plot4), file(plot5), file(plot6), file(plot7), file(plot8), file(plot9), file(plot10), file(plot11), file(plot12), file(plot13), file(plot14) 

    output:
    tuple val(meta), path("*_allSNVdiagnosticsPlots*.pdf") , emit: plots
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rerun  = params.rerunfiltering ? "1": ''

    """
    export TEMP=\$(mktemp -d)
    export TMPDIR=\$TEMP
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \\
        -sOutputFile=${prefix}_allSNVdiagnosticsPlots${rerun}.pdf \\
        $plot1 $plot2 $plot3 $plot4 $plot5 $plot6 $plot7 $plot8 $plot9 $plot10 $plot11 $plot12 $plot13 $plot14
    rm -rf \$TEMP
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ghostscript: \$(echo \$(gs -v 2>&1) | sed 's/^.*GPL Ghostscript //; s/ .*\$//')
    END_VERSIONS
    """ 
}
