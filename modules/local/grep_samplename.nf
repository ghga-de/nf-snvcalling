process GREP_SAMPLENAME {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/samtools:v1.9':'kubran/samtools:v1.9' }"

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
        controlname=`samtools view -H $control | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`
        tumorname=`samtools view -H $tumor | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS

        """
    }
    else {
        """
        controlname='dummy'
        tumorname=`samtools view -H $tumor | grep '^@RG' | sed "s/.*SM:\\([^\\t]*\\).*/\\1/g" | uniq`
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """
    }
}
