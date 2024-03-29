//Annotate with polymorphisms (dbSNP, 1K genomes, and local controls) and prepare annovar input file
process ANNOTATE_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple path(kgenome), path(kgenome_i) 
    tuple path(dbsnpsnv), path(dbsnpsnv_i) 
    tuple path(localcontrolwgs), path(localcontrolwgs_i) 
    tuple path(localcontrolwes), path(localcontrolwes_i) 
    tuple path(gnomadgenomes), path(gnomadgenomes_i) 
    tuple path(gnomadexomes), path(gnomadexomes_i)
    val (chrprefix)

    output:
    tuple val(meta), path('*.ForAnnovar.bed')    , emit: forannovar
    tuple val(meta), path('*.vcf')               , emit: unziped_vcf
    tuple val(meta), path("*_median.txt")        , emit: median
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def cmdfilter   = meta.iscontrol == "1" ? "| median.pl - vcf_control_median.txt" : ""
    def pipe  = ["${cmdfilter}",
                dbsnpsnv ? " | annotate_vcf.pl -a - -b ${dbsnpsnv} --columnName='DBSNP' --reportMatchType --bAdditionalColumn=2  --reportLevel 4" : '',
                kgenome ? " | annotate_vcf.pl -a - -b ${kgenome} --columnName='1K_GENOMES' --reportMatchType --bAdditionalColumn=2 --reportLevel 4" : '',
                localcontrolwgs ? " | annotate_vcf.pl -a - -b ${localcontrolwgs} --columnName='LocalControlAF_WGS' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                localcontrolwes ? " | annotate_vcf.pl -a - -b ${localcontrolwes} --columnName='LocalControlAF_WES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadgenomes ? " | annotate_vcf.pl -a - -b ${gnomadgenomes} --columnName='GNOMAD_GENOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : '',
                gnomadexomes ? " | annotate_vcf.pl -a - -b ${gnomadexomes} --columnName='GNOMAD_EXOMES' --bFileType vcf --reportLevel 4 --reportMatchType" : ''
                ].join(' ').trim()

    def maxcontrolcov = meta.iscontrol == "1" ? "[[ `cat vcf_control_median.txt` -lt 1 ]] && echo 'Median was not calculated correctly' && exit 3": "touch vcf_nocontrol_median.txt"
    """
    zcat < $vcf $pipe | \\
        tee ${prefix}.vcf | vcf_to_annovar.pl $chrprefix "" > ${prefix}.ForAnnovar.bed

    $maxcontrolcov 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}