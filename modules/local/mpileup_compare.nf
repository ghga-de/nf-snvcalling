process MPILEUP_COMPARE {
    tag "$meta.id $intervals"
    label 'process_single'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/samtools:v1.9':'kubran/samtools:v1.9' }"

    input:
    tuple val(meta), path(vcf), val(intervals)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.nppileup.vcf")     , emit: vcf
    path  "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_c = intervals == "contig" ? "$args" : "$args -r ${intervals}"

    if (meta.iscontrol == '1' && params.runCompareGermline)
    {
        """
        samtools mpileup \\
            $args_c \\
            -l $vcf \\
            -f $fasta \\
            $meta.control_bam > ${prefix}.${intervals}.control.temp
        
        vcf_pileup_compare_allin1_basecount.pl $vcf \\
            ${prefix}.${intervals}.control.temp "Header" > ${prefix}.${intervals}.nppileup.vcf        

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """

    }
    else {
        """
        mv $vcf ${prefix}.${intervals}.nocontrol.nppileup.vcf 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """
    }
}
