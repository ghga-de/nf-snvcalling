process FILTER_STRAND_BIAS {
    tag "$meta.id"
    label 'process_intermediate'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'kubran/odcf_snvcalling:v0':'odcf_snvcalling_v0.sif' }"

    publishDir params.outdir+'/mpileup' , mode: 'copy'

    debug true

    input:
    tuple val(meta), val(intervals), path(vcf)
    tuple path(fasta), path(fai) 

    output:
    tuple val(meta), val(intervals), path("*.bias.vcf")     , emit: vcf
    path  "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    firstline=`cat $vcf | grep -v "^#" | head -n1`
    if [[ -n "$firstLineVCF" ]]; then
        seqConText_annotator.pl $vcf $fasta 10 | \\
        filterVcfForBias.py --outf=${prefix}.${intervals}.bias.vcf $args

    else
        mv $vcf ${prefix}.${intervals}.empty.bias.vcf
        echo "WARNING! Found no variant at ${intervals} "
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
