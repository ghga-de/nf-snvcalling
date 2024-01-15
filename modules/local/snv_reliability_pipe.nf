
process SNV_RELIABILITY_PIPE {
    tag "$meta.id"
    label 'process_single'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta),path(ch_vcf),path(ch_vcf_i)
    tuple path(repeatmasker),path(repeatmasker_i)
    tuple path(dacblacklist),path(dacblacklist_i)
    tuple path(dukeexcluded),path(dukeexcluded_i)
    tuple path(hiseqdepth),path(hiseqdepth_i)
    tuple path(selfchain),path(selfchain_i)
    tuple path(mapability),path(mapability_i)
    tuple path(simpletandemrepeats),path(simpletandemrepeats_i)

    output:
    tuple val(meta), path('*.annotated.vcf')   , emit: vcf
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def pipe       = [repeatmasker ? " | annotate_vcf.pl -a - -b ${repeatmasker} --bFileType=bed --columnName='REPEAT_MASKER'" : '',
                    dacblacklist ? " | annotate_vcf.pl -a - -b ${dacblacklist}  --bFileType=bed --columnName='DAC_BLACKLIST'" : '',
                    dukeexcluded ? " | annotate_vcf.pl -a - -b ${dukeexcluded} --bFileType=bed --columnName='DUKE_EXCLUDED'" : '',
                    hiseqdepth ? " | annotate_vcf.pl -a - -b ${hiseqdepth} --bFileType=bed --columnName='HISEQDEPTH'" : '',
                    selfchain ? " | annotate_vcf.pl -a - -b ${selfchain} --bFileType=bed --columnName='SELFCHAIN' --bAdditionalColumns=4 --maxNrOfMatches=5" : '',
                    mapability ? " | annotate_vcf.pl -a - -b ${mapability} --bFileType=bed --columnName='MAPABILITY' --breakPointMode --aEndOffset=1" : '',
                    simpletandemrepeats ? " | annotate_vcf.pl -a - -b ${simpletandemrepeats} --bFileType=bed --columnName='SIMPLE_TANDEMREPEATS' --bAdditionalColumns=4" : ''
                    ].join(' ').trim()

    """
    zcat < $ch_vcf $pipe > ${prefix}.annotated.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}