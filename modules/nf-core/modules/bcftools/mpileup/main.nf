process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai), val(tumorname), val(controlname)
    tuple path(fasta), path(fai) 
    val(intervals)


    output:
    tuple val(meta), path("*.gz")      , emit: vcf
    tuple val(meta), path("*.tbi")     , emit: tbi
    tuple val(meta), path("*stats.txt"), emit: stats
    tuple val(meta), path("*.mpileup") , emit: mpileup, optional: true
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        $args \\
        $tumor \\
        -r $intervals \\
        | bcftools call --output-type v $args2 \\
        | bcftools view --output-file ${prefix}.vcf.gz --output-type z $args3

    tabix -p vcf -f ${prefix}.vcf.gz

    bcftools stats ${prefix}.vcf.gz > ${prefix}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}