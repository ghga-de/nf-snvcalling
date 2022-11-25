//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATION_PIPES {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'library://kubran/odcf/odcf_platypusindelcalling:v0' :'kubran/odcf_platypusindelcalling:v0' }"

    debug true
    
    input:
    tuple val(meta)           , file(vcf)             , file(vcf_tbi)
    tuple file(enchangers)    , file(enchangers_i)
    tuple file(cpgislands)    , file(cpgislands_i)
    tuple file(tfbscons)      , file(tfbscons_i)
    tuple file(encode_dnase)  , file(encode_dnase_i)
    tuple file(mirnas_snornas), file(mirnas_snornas_i)
    tuple file(cosmic)        , file(cosmic_i)
    tuple file(mirbase)       , file(mirbase_i)
    tuple file(mir_targets)   , file(mir_targets_i)
    tuple file(cgi_mountains) , file(cgi_mountains_i)
    tuple file(phastconselem) , file(phastconselem_i)
    tuple file(encode_tfbs)   , file(encode_tfbs_i)

    output:
    tuple val(meta), path('*.deepanno.vcf.gz'), path('*.deepanno.vcf.gz.tbi') , emit: vcf
    path  "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def pipe  = [enchangers.baseName !='input' ? " | annotate_vcf.pl -a - -b ${enchangers} --bFileType=bed --columnName='Enhancers'" : '',
                cpgislands.baseName !='input' ? " | annotate_vcf.pl -a - -b ${cpgislands} --bFileType=bed --columnName='CpGislands'" : '',
                tfbscons.baseName !='input' ? " | annotate_vcf.pl -a - -b ${tfbscons} --bFileType=bed --columnName='TFBScons'" : '',
                mirnas_snornas.baseName !='input' ? " | annotate_vcf.pl -a - -b ${mirnas_snornas} --bFileType=bed --columnName='miRNAs_snoRNAs'" : '',
                encode_dnase.baseName !='input' ? " | annotate_vcf.pl -a - -b ${encode_dnase} --bFileType=bed --columnName='ENCODE_DNASE'" : '',
                mirbase.baseName !='input' ? " | annotate_vcf.pl -a - -b ${mirbase} --bFileType=bed --columnName='miRBase18'" : '',
                cosmic.baseName !='input' ? " | annotate_vcf.pl -a - -b ${cosmic} --bFileType=bed --columnName='COSMIC' --bAdditionalColumns=7,8,9 --reportLevel=1" : '',
                mir_targets.baseName !='input' ? " | annotate_vcf.pl -a - -b ${mir_targets} --columnName='miRNAtargets'" : '' ,
                cgi_mountains.baseName !='input' ? " | annotate_vcf.pl -a - -b ${cgi_mountains} --bFileType=bed --columnName='CgiMountains' --bAdditionalColumns=4" : '',
                phastconselem.baseName !='input' ? " | annotate_vcf.pl -a - -b ${phastconselem} --bFileType=bed --columnName='phastConsElem20bp' --bAdditionalColumns=4" : '',
                encode_tfbs.baseName !='input' ? " | annotate_vcf.pl -a - -b ${encode_tfbs} --columnName='ENCODE_TFBS'" : ''
                ].join(' ').trim() 
    """
    zcat < $vcf $pipe > ${prefix}.deepanno.vcf
    bgzip -c ${prefix}.deepanno.vcf > ${prefix}.deepanno.vcf.gz
    tabix ${prefix}.deepanno.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}