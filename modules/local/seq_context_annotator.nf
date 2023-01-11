process SEQ_CONTEXT_ANNOTATOR {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"

    input:
    tuple val(meta), path(vcf), val(intervals)
    tuple path(fasta), path(fai) 

    output:
    tuple val(meta), path("*.bias.vcf")                , emit: vcf
    tuple val(meta), path("*.bias.bcftools_stats.txt") , emit: stats 
    tuple val(meta), val(intervals)                    , emit: intervals 
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqContext_annotator.pl fastaFromBed $vcf $fasta 10 | \\
        rawSnvFilter.py --outf=${prefix}.${intervals}.bias.vcf $args

    bcftools stats ${prefix}.${intervals}.bias.vcf > ${prefix}.${intervals}.bias.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
