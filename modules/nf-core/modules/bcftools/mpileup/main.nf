process BCFTOOLS_MPILEUP {
    tag "$meta.id $intervals"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bcftools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.9--h47928c2_2':
        'quay.io/biocontainers/bcftools:1.9--h47928c2_2' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai), val(tumorname), val(controlname), val(intervals), path(interval_file)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf")               , emit: vcf
    tuple val(meta), path("*.bcftools_stats.txt"), emit: stats 
    tuple val(meta), val(intervals)          , emit: intervals 
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def args2    = task.ext.args2 ?: ''
    def args3    = task.ext.args3 ?: ''
    def prefix   = task.ext.prefix ?: "${meta.id}"
    def args_c   = interval_file ? "$args2 -R ${interval_file}" : "$args -r ${intervals}"
    def ref_spec = params.fasta.contains("38") ? "$args3 --ploidy GRCh38": "$args3"

    """
    bcftools \\
        mpileup \\
        --fasta-ref $fasta \\
        $args_c \\
        $tumor \\
        | bcftools call --output-type v $ref_spec > ${prefix}.${intervals}.vcf

    bcftools stats ${prefix}.${intervals}.vcf > ${prefix}.${intervals}.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
