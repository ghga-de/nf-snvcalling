//purityEST - from Florian. Needs the original SNV file because it also considers germline (DP5 field)
//has everything hardcoded (in which fields to look and confidence 8)
process PURITY_RELOADED {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), file(vcf), file(index)

    output:
    tuple val(meta), path('*.txt')        , emit: purity
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    PurityReloaded.py \\
        $vcf \\
        `findConfColumn.pl ${vcf}` > ${prefix}_purityEST.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}