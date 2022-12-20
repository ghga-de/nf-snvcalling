//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"

    input:
    tuple val(meta)            , file(vcf)     , file(vcf_tbi), val(tumorname), val(controlname)
    tuple file(kgenome)        , file(kgenome_i)
    tuple file(dbsnpindel)     , file(dbsnpindel_i)
    tuple file(exac)           , file(exac_i)
    tuple file(evs)            , file(evs_i)
    tuple file(localcontrolwgs), file(localcontrolwgs_i)
    tuple file(localcontrolwes), file(localcontrolwes_i)
    tuple file(gnomadgenomes)  , file(gnomadgenomes_i)
    tuple file(gnomadexomes)   , file(gnomadexomes_i)
    val (chrprefix)

    output:
    tuple val(meta), path('*.ForAnnovar.bed')                         , emit: forannovar
    tuple val(meta), path('*.vcf')                                    , emit: unziped_vcf
    path("*.vcf_control_median.txt")                                  , optional: true
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def cmdfilter   = meta.iscontrol == "1" ? "| median.pl - ${prefix}.vcf_control_median.txt" : ""
    def pipe  = ["${cmdfilter}",
                kgenome.baseName !='input' ? " | annotate_vcf.pl -a - -b ${dbsnpindel} --columnName='DBSNP' --reportMatchType --bAdditionalColumn=2  --reportLevel 4" : '',
                dbsnpindel.baseName !='input' ? " | annotate_vcf.pl -a - -b ${kgenome} --columnName='1K_GENOMES' --reportMatchType --bAdditionalColumn=2 --reportLevel 4" : '',
                exac.baseName !='input' ? " | annotate_vcf.pl -a - -b ${exac} --columnName='ExAC' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                evs.baseName !='input' ? " | annotate_vcf.pl -a - -b ${evs} --columnName='EVS' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                localcontrolwgs.baseName !='input' ? " | annotate_vcf.pl -a - -b ${localcontrolwgs} --columnName='LocalControlAF_WGS' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                localcontrolwes.baseName !='input' ? " | annotate_vcf.pl -a - -b ${localcontrolwes} --columnName='LocalControlAF_WES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadgenomes.baseName !='input' ? " | annotate_vcf.pl -a - -b ${gnomadgenomes} --columnName='GNOMAD_GENOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadexomes.baseName !='input' ? " | annotate_vcf.pl -a - -b ${gnomadexomes} --columnName='GNOMAD_EXOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : ''
                ].join(' ').trim()

    """
    zcat < $vcf $pipe |\\
        tee ${prefix}.vcf | vcf_to_annovar.pl $chrprefix "" > ${prefix}.ForAnnovar.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}