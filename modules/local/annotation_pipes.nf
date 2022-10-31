//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATION_PIPES {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'odcf_indelcalling_v5.sif' :
    'kubran/odcf_indelcalling:v5' }"


    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    path 'deneme.txt'
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def enchangers_seq     = params.enchangers ? "-en ${params.enchangers}" : ''
    def cpgislands_seq     = params.cpgislands ? "-cp ${params.cpgislands}" : ''
    def tfbscons_seq       = params.tfbscons ? "-tf ${params.tfbscons}" : ''
    def mirnas_snornas_seq = params.mirnas_snornas ? "-ms ${params.mirnas_snornas}" : ''


    def pipe               = [params.enchangers ? "-en ${params.enchangers}" : '',
                            params.cpgislands ? "-cp ${params.cpgislands}" : '',
                            params.tfbscons ? "-tf ${params.tfbscons}" : '',
                            params.mirnas_snornas ? "-ms ${params.mirnas_snornas}" : ''].join(' ').trim()

    """
    cat $tumor $control $pipe > deneme.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
