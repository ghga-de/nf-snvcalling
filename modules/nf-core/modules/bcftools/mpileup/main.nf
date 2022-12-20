process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1':
        'quay.io/biocontainers/bcftools:1.16--hfe4b78e_1' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai), val(tumorname), val(controlname), val(intervals)
    tuple path(fasta), path(fai) 

    output:
    tuple val(meta), val(intervals), path("*.vcf")               , emit: vcf
    tuple val(meta), val(intervals), path("*.bcftools_stats.txt"), emit: stat 
    path  "versions.yml"                                         , emit: versions

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
        | bcftools call --output-type v $args2 > ${prefix}.${intervals}.vcf

    bcftools stats ${prefix}.${intervals}.vcf > ${prefix}.${intervals}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
