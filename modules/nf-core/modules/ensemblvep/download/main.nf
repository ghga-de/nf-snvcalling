process ENSEMBLVEP_DOWNLOAD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.0--pl5321h2a3209d_0' :
        'quay.io/biocontainers/ensembl-vep:110.0--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(x)

    output:
    path("vep_cache")      , emit: cache
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    vep_install \\
        --CACHEDIR vep_cache \\
        --SPECIES "homo_sapiens" \\
        --ASSEMBLY $params.vep_genome \\
        --CACHE_VERSION $params.vep_cache_version \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir vep_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
