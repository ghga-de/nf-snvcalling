process SORT_NONSTANDARD_VCF {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/samtools:v1.9':'kubran/samtools:v1.9' }"

    input:
    tuple val(meta), path(vcf_gz), path(index)

    output:
    tuple val(meta), path("*_raw.vcf.gz"),path("*_raw.vcf.gz.tbi")  , emit: output
    path "versions.yml"                                             , emit: versions

    script: 
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    (zcat $vcf_gz | head -n 5000 | \\
        grep "#" ; zcat $vcf_gz | \\
        grep -v "#" | \\
        sort -T . -k1,1V -k2,2n ) | bgzip > snvs_${prefix}_raw.vcf.gz

    tabix -p vcf snvs_${prefix}_raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
    END_VERSIONS
    """
    
}