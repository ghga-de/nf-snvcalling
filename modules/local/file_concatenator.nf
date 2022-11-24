process FILE_CONCATENATOR {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'kubran/odcf_snvcalling:v2':'odcf_snvcalling_v2.sif' }"

    input:
    tuple val(meta), path(vcfs)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi")      , emit: vcf
    path  "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    headeredFileConcatenator.pl $vcfs > snv_${prefix}.vcf

    bgzip snv_${prefix}.vcf && tabix -p vcf snv_${prefix}.vcf.gz 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS

    """ 
}
