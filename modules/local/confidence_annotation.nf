//MAX_CONTROL_COV=0 is not performed!
process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'library://kubran/odcf/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple val(meta), file(a), file(b), val(tumorname), val(controlname)

    output:
    tuple val(meta), path('*.vcf.gz') ,   path('*.vcf.gz.tbi')    , emit: vcf_ann
    path  "versions.yml"                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag    = meta.iscontrol == "1" ? "" : "--nocontrol"
    def samples    = meta.iscontrol == "1" ? "--controlColName=$controlname --tumorColName=$tumorname" : "--nocontrol --tumorColName=$tumorname"

    """
    cat > $vcfgz | confidenceAnnotation_SNVs.py $controlflag -i - $args -a 0 > snv_${prefix}.conf.vcf

    bgzip -c snv_${prefix}.conf.vcf > snv_${prefix}.conf.vcf.gz
    tabix snv_${prefix}.conf.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS

    """
    
}