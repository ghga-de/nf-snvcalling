process MPILEUP_COMPARE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'kubran/odcf_snvcalling:v2':'odcf_snvcalling_v2.sif' }"

    debug true

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

    if (meta.iscontrol == '1' && params.runCompareGermline)
    {
        """
        samtools mpileup $args -r $intervals -l $vcf -f $fasta $meta.control_bam > ${prefix}.${intervals}.control.temp
        
        vcf_pileup_compare_allin1_basecount.pl $vcf ${prefix}.${intervals}.control.temp "Header" > ${prefix}.${intervals}.nppileup.vcf        

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS

        """

    }
    else {
        """
        mv $vcf ${prefix}.${intervals}.nocontrol.nppileup.vcf 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
