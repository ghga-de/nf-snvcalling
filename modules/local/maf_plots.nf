process MAF_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(somatic_snv), file(maf_values), file(indbsnp)

    output:
    tuple val(meta), path('*.pdf')       , emit: plot
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    MAF_plots.r \\
        $maf_values \\
        `grep -v "^#" ${somatic_snv} | wc -l` \\
        ${prefix}_MAF_conf_${params.min_confidence_score}_to_10.pdf \\
        $prefix \\
        ` awk '{FS="\t"}{if(NR==2)print \$5}'	${indbsnp}`

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}