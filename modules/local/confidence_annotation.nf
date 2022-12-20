//MAX_CONTROL_COV=0 is not performed!
process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)

    output:
    tuple val(meta),path('*.vcf.gz') ,path('*.vcf.gz.tbi')    , emit: vcf_ann
    path  "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag    = meta.iscontrol == "1" ? "" : "--nocontrol"
    def ref_hg38   = params.ref_type == "hg38" ? "--refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "" 
    def ref_spec   = params.ref_type == "hg38" ? "--gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf}" : ""
    """
    cat > $vcfgz | confidenceAnnotation_SNVs.py $controlflag -i - $args $ref_hg38 -a 0 $ref_spec > snv_${prefix}.conf.vcf

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