//snv_reliability_pipe

process SNV_RELIABILITY_PIPE {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'library://kubran/odcf/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"
    
    input:
    tuple val(meta),                 file(ch_vcf),               file(ch_vcf_i)
    tuple path(repeatmasker),        path(repeatmasker_i)
    tuple path(dacblacklist),        path(dacblacklist_i)
    tuple path(dukeexcluded),        path(dukeexcluded_i)
    tuple path(hiseqdepth),          path(hiseqdepth_i)
    tuple path(selfchain),           path(selfchain_i)
    tuple path(mapability),          path(mapability_i)
    tuple path(simpletandemrepeats), path(simpletandemrepeats_i)

    output:
    tuple val(meta), path('*.annotated.vcf.gz'), path('*.annotated.vcf.gz.tbi')   , emit: vcf
    path  "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"

    """
    zcat < $ch_vcf | \\
        annotate_vcf.pl -a - -b $repeatmasker --bFileType=bed --columnName='REPEAT_MASKER' | \\
        annotate_vcf.pl -a - -b $dacblacklist --bFileType=bed --columnName='DAC_BLACKLIST' | \\
        annotate_vcf.pl -a - -b $dukeexcluded --bFileType=bed --columnName='DUKE_EXCLUDED' | \\
        annotate_vcf.pl -a - -b $hiseqdepth --bFileType=bed --columnName='HISEQDEPTH' | \\
        annotate_vcf.pl -a - -b $selfchain --bFileType=bed --columnName='SELFCHAIN' --bAdditionalColumns=4 --maxNrOfMatches=5 | \\
        annotate_vcf.pl -a - -b $mapability --bFileType=bed --columnName='MAPABILITY' | \\
        annotate_vcf.pl -a - -b $simpletandemrepeats --bFileType=bed --columnName='SIMPLE_TANDEMREPEATS' --bAdditionalColumns=4 \\
        > ${prefix}.annotated.vcf

    bgzip -c ${prefix}.annotated.vcf > ${prefix}.annotated.vcf.gz
    tabix ${prefix}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}