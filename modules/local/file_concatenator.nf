process MPILEUP_COMPARE {
    tag "$meta.id"
    label 'process_intermediate'

    conda (params.enable_conda ? "bioconda::bcftools=1.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'kubran/odcf_snvcalling:v0':'odcf_snvcalling_v0.sif' }"

    publishDir params.outdir+'/mpileup' , mode: 'copy'

    debug true

    input:
    tuple val(meta),  val(intervals), path(vcf)
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai), val(tumorname), val(controlname), val(intervals)
    tuple path(fasta), path(fai) 

    output:
    tuple val(meta), path("*.tmp")     , emit: vcf
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.iscontrol == '1' && params.runCompareGermline)
    {
  
        """
        firstline=`cat $vcf | grep -v "^#" | head -n1`

        if [[ -n "$firstLineVCF" ]]; then
            samtools mpileup $args -r $intervals -l $vcf -f $fasta $control > \\              
                NP_MPILEUP${intervals} & \\
                vcf_pileup_compare_allin1_basecount.pl $vcf NP_MPILEUP${intervals} "Header" > \\
                ${prefix}.nppileup.${intervals}.vcf
        else
            mv $vcf ${prefix}.${intervals}.empty.vcf
        fi          

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS

        """

    }
    else {
        """
        controlname='dummy'
        tumorname=`samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }
}
