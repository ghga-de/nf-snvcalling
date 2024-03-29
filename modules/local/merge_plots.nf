process MERGE_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(plot1), path(plot2)

    output:
    tuple val(meta), path("*_allSNVdiagnosticsPlots*.pdf") , emit: plots
    path  "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    export TEMP=\$(mktemp -d)
    export TMPDIR=\$TEMP
    gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite \\
        -sOutputFile=snvs_${prefix}_allSNVdiagnosticsPlots.pdf \\
        $plot1 $plot2
    rm -rf \$TEMP
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ghostscript: \$(echo \$(gs -v 2>&1) | sed 's/^.*GPL Ghostscript //; s/ .*\$//')
    END_VERSIONS
    """ 
}
