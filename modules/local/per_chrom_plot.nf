process PER_CHROM_PLOT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"
    
    input:
    tuple val(meta), file(distance)
    each file(chr_file)

    output:
    tuple val(meta), path('*.pdf')   , emit: plot        
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def rerun      = params.rerunfiltering ? "_filteredAltMedian${params.median_filter_threshold}": ""     
    
    """
    snvsPerChromPlot.r \\
        -i $distance \\
        -l $chr_file \\
        -s ${prefix} \\
        -o ${prefix}_perChromFreq_conf_${params.min_confidence_score}_to_10${rerun}.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}