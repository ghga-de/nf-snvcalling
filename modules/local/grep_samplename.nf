process GREP_SAMPLENAME {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v7.sif' :'kubran/odcf_indelcalling:v7' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    tuple val(meta), env(tumorname)    , env(controlname)          , emit: samplenames
    path "versions.yml"     , emit: versions

    script: 
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    if (meta.iscontrol == '1')
    {
        """
        controlname=`samtools view -H $meta.control_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`
        tumorname=`samtools view -H $meta.tumor_bam | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`

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
