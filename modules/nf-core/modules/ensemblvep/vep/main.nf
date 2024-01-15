process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:110.0--pl5321h2a3209d_0' :
        'quay.io/biocontainers/ensembl-vep:110.0--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(vcf)
    path(cache)
    tuple path(fasta), path(index)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi")  , emit: vcf
    path "*.summary.html"                                    , emit: report
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta $fasta" : ""
    def dir_cache = args.contains("--offline") ? "--dir_cache ${cache} --cache" : "--database"
    """
    vep \\
        -i $vcf \\
        -o ${prefix}.vep.vcf.gz \\
        --format vcf \\
        --vcf \\
        --compress_output bgzip \\
        $args \\
        $reference \\
        --assembly $params.vep_genome \\
        --species "homo_sapiens" \\
        --cache_version $params.vep_cache_version \\
        $dir_cache \\
        --stats_file ${prefix}.summary.html

        tabix ${prefix}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.tab.gz
    touch ${prefix}.json.gz
    touch ${prefix}.summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}
