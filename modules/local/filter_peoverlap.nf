//# create BaseScore FIFOs and their consumer processes (zip and write to target file)
//# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}

process FILTER_PEOVERLAP {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*.allelebasescore.vcf.gz'), path('*.allelebasescore.vcf.gz.tbi')   , emit: allelebasescore_vcf
    tuple val(meta), path('*somatic_snvs_for_bias.vcf')                                       , emit: somatic_snvs_for_bias   
    path  "versions.yml"                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def args2      = task.ext.args2 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag = meta.iscontrol == "1" ? "" : "--nocontrol"
    def ref_hg38   = params.ref_type == "hg38" ? "--refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "" 
    def ref_spec   = params.ref_type == "hg38" ? "--gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf}" : ""

    if (params.runArtifactFilter){
        """
        cat > $vcfgz | filter_PEoverlap.py \\
            $controlflag \\
            --alignmentFile=$meta.tumor_bam \\
            $args \\
            --altBaseQualFile=${prefix}_alternative_allele_base_qualities.txt \\
            --refBaseQualFile=${prefix}_reference_allele_base_qualities.txt \\
            --altBasePositionsFile=${prefix}_alternative_allele_read_positions.txt \\
            --refBasePositionsFile=${prefix}_reference_allele_read_positions.txt \\
            --referenceFile=$fasta | \\
                confidenceAnnotation_SNVs.py $controlflag -i - $args2 $ref_hg38 -a 0 -f ${prefix}_somatic_snvs_for_bias.vcf \\
                    $ref_spec > snv_${prefix}.allelebasescore.vcf

        bgzip -c snv_${prefix}.allelebasescore.vcf > snv_${prefix}.allelebasescore.vcf.gz
        tabix snv_${prefix}.allelebasescore.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
        END_VERSIONS
        """   
    }
    else {
        """
        cat > $vcfgz | filter_PEoverlap.py \\
            $controlflag \\
            --alignmentFile=$meta.tumor_bam \\
            $args \\
            --altBaseQualFile=${prefix}_alternative_allele_base_qualities.txt \\
            --refBaseQualFile=${prefix}_reference_allele_base_qualities.txt \\
            --altBasePositionsFile=${prefix}_alternative_allele_read_positions.txt \\
            --refBasePositionsFile=${prefix}_reference_allele_read_positions.txt \\
            --referenceFile=$fasta | \\
                confidenceAnnotation_SNVs.py $controlflag -i - $args2 $ref_hg38 \\
                    $ref_spec > snv_${prefix}.allelebasescore.vcf

        bgzip -c snv_${prefix}.allelebasescore.vcf > snv_${prefix}.allelebasescore.vcf.gz
        tabix snv_${prefix}.allelebasescore.vcf.gz

        touch ${prefix}_somatic_snvs_for_bias.vcf 

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
        END_VERSIONS
        """

    }
}