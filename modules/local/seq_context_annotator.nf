process SEQ_CONTEXT_ANNOTATOR {
    tag "$meta.id $intervals"
    label 'process_single'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/samtools:v1.9':'kubran/samtools:v1.9' }"

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
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    seqContext_annotator.pl fastaFromBed $vcf $fasta 10 | \\
        rawSnvFilter.py --outf=${prefix}.${intervals}.bias.vcf $args

    bcftools stats ${prefix}.${intervals}.bias.vcf > ${prefix}.${intervals}.bias.bcftools_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
