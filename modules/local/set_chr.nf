process SET_CHR {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v7.sif' :'kubran/odcf_indelcalling:v7' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    env(chr_prefix)                 , emit: chr
    path "versions.yml"             , emit: versions

    script: 
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    #!/bin/bash
    set -o pipefail

    x=`samtools view -H $tumor | grep "^@SQ" | grep "SN:chr" | wc -l`

    if [[ x -gt 0 ]]; then
        chr_prefix="chr"
    else
        chr_prefix="dummy"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}