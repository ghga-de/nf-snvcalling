process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::bcftools=1.16" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'kubran/odcf_snvcalling:v0':'odcf_snvcalling_v0.sif'}"

    publishDir params.outdir+'/mpileup' , mode: 'copy'

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai), val(tumorname), val(controlname), val(intervals)
    tuple path(fasta), path(fai) 

    output:
    tuple val(meta), val(intervals), path("*.vcf")     , emit: vcf
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools mpileup --fasta-ref $fasta --regions $intervals $args $tumor | bcftools call $args2 > ${prefix}.${intervals}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
