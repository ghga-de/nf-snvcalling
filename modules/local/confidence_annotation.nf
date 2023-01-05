process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), file(vcf) file(median)

    output:
    tuple val(meta),path('*.conf.vcf.gz') ,path('*.conf.vcf.gz.tbi')    , emit: vcf_ann
    path  "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag    = meta.iscontrol == "1" ? "" : "--nocontrol"
    def confoptions = meta.iscontrol == "1" ? "-t `cat $medium`" : "$args"
    def ref_hg38   = params.ref_type == "hg38" ? "--refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "" 
    def ref_spec   = params.ref_type == "hg38" ? "--gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf}" : ""

    """
    mkfifo snvAnnotationFIFO_${prefix}.vcf
    if [[ "$meta.iscontrol" == "1" ]]; then
        cat $vcf > ${} &
    else
        cat $vcf | confidenceAnnotation_SNVs.py $controlflag -i - $confoptions  $ref_hg38 -a 0 $ref_spec > snvAnnotationFIFO_${prefix}.vcf &
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
    
}