//Annotate with polymorphisms (dbSNP, 1K genomes, ExAC, EVS, and local controls) and prepare annovar input file
process ANNOTATION_PIPES {
    tag "$meta.id"
    label 'process_single'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(vcf), path(vcf_tbi)
    tuple path(enchangers),path(enchangers_i)
    tuple path(cpgislands),path(cpgislands_i)
    tuple path(tfbscons),path(tfbscons_i)
    tuple path(encode_dnase),path(encode_dnase_i)
    tuple path(mirnas_snornas),path(mirnas_snornas_i)
    tuple path(cosmic),file(cosmic_i)
    tuple path(mirbase),path(mirbase_i)
    tuple path(mir_targets),path(mir_targets_i)
    tuple path(cgi_mountains),path(cgi_mountains_i)
    tuple path(phastconselem),path(phastconselem_i)
    tuple path(encode_tfbs),path(encode_tfbs_i)
    tuple path(mirnas_sncrnas),path(mirnas_sncrnas_i)

    output:
    tuple val(meta), path('*.vcf.gz'), path('*.vcf.gz.tbi') , emit: vcf
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def pipe  = [enchangers ? " | annotate_vcf.pl -a - -b ${enchangers} --bFileType=bed --columnName='Enhancers'" : '',
                cpgislands ? " | annotate_vcf.pl -a - -b ${cpgislands} --bFileType=bed --columnName='CpGislands'" : '',
                tfbscons ? " | annotate_vcf.pl -a - -b ${tfbscons} --bFileType=bed --columnName='TFBScons'" : '',
                mirnas_snornas ? " | annotate_vcf.pl -a - -b ${mirnas_snornas} --bFileType=bed --columnName='miRNAs_snoRNAs'" : '',
                encode_dnase ? " | annotate_vcf.pl -a - -b ${encode_dnase} --bFileType=bed --columnName='ENCODE_DNASE'" : '',
                mirbase ? " | annotate_vcf.pl -a - -b ${mirbase} --bFileType=bed --columnName='miRBase'" : '',
                cosmic ? " | annotate_vcf.pl -a - -b ${cosmic} --bFileType=bed --columnName='COSMIC' --bAdditionalColumns=7,8,9 --reportLevel=1" : '',
                mir_targets ? " | annotate_vcf.pl -a - -b ${mir_targets} --columnName='miRNAtargets'" : '' ,
                cgi_mountains ? " | annotate_vcf.pl -a - -b ${cgi_mountains} --bFileType=bed --columnName='CgiMountains' --bAdditionalColumns=4" : '',
                phastconselem ? " | annotate_vcf.pl -a - -b ${phastconselem} --bFileType=bed --columnName='phastConsElem20bp' --bAdditionalColumns=4" : '',
                encode_tfbs ? " | annotate_vcf.pl -a - -b ${encode_tfbs} --columnName='ENCODE_TFBS'" : '',
                mirnas_sncrnas ? " | annotate_vcf.pl -a - -b ${mirnas_sncrnas} --bFileType=bed --columnName='miRNAs_sncRNAs'" : ''
                ].join(' ').trim() 
    """
    zcat < $vcf $pipe > snvs_${prefix}.vcf
    bgzip -c snvs_${prefix}.vcf > snvs_${prefix}.vcf.gz
    tabix snvs_${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}