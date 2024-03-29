process CONFIDENCE_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"
    
    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path('*.confidence.vcf')    , emit: vcf
    path  "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def confoptions = params.fasta.contains("38") ? "${params.confidenceoptions} --refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "${params.confidenceoptions}" 

    if (meta.iscontrol == "1"){
        """
        mv $vcf snv_${prefix}.confidence.vcf
        cat snv_${prefix}.confidence.vcf > ${}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python2 --version 2>&1 | sed 's/Python //g')
        END_VERSIONS
    """
    } 
    else{
        """
        cat < $vcf | confidenceAnnotation_SNVs.py \\
            --nocontrol \\
            -i - \\
            -a 0 \\
            -p \\
            $confoptions \\
            --gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} \\
            --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} \\
            --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf} > snv_${prefix}.confidence.vcf       

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python2 --version 2>&1 | sed 's/Python //g')
        END_VERSIONS
        """
    }
}