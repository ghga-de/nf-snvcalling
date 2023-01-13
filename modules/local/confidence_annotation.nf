process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(vcf)

    output:
    tuple val(meta), path('*.confidence.vcf')    , emit: vcf
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def confoptions   = params.ref_type == "hg38" ? "${params.confidenceoptions} --refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "${params.confidenceoptions}" 

    if (meta.iscontrol == "1"){
        """
        mv $vcf snv_${prefix}.confidence.vcf
        cat snv_${prefix}.confidence.vcf > ${}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
    """
    } 
    else{
        """
        cat < $vcf | confidenceAnnotation_SNVs.py \\
            --nocontrol \\
            -i - \\
            -a 0 \\
            $confoptions \\
            --gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} \\
            --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} \\
            --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf} > snv_${prefix}.confidence.vcf       

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: v5.28.1
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
}