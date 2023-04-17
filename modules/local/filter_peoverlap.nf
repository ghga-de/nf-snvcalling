//# create BaseScore FIFOs and their consumer processes (zip and write to target file)
//# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}
process FILTER_PEOVERLAP {
    tag "$meta.id"
    label 'process_medium_high'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/samtools:v1.9':'kubran/samtools:v1.9' }"
    
    input:
    tuple val(meta),   file(vcf)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*_peoverlap.vcf')                         , emit: vcf
    tuple val(meta), path('*_somatic_snvs_for_bias.vcf')             , emit: somatic_snvs
    tuple val(meta), path('*_alternative_allele_base_qualities.txt') , emit: alternative_allele_base_qualities
    tuple val(meta), path('*_reference_allele_base_qualities.txt')   , emit: reference_allele_base_qualities 
    tuple val(meta), path('*_alternative_allele_read_positions.txt') , emit: alternative_allele_read_positions 
    tuple val(meta), path('*_reference_allele_read_positions.txt')   , emit: reference_allele_read_positions
    path  "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def args2      = task.ext.args2 ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag   = meta.iscontrol == "1" ? "" : "--nocontrol "
    def confoptions   = params.fasta.contains("38") ? "${params.confidenceoptions} --refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "${params.confidenceoptions}"

    """
    cat < $vcf | filter_PEoverlap.py $controlflag\\
        --alignmentFile=$meta.tumor_bam \\
        --mapq=${params.mapqual} \\
        --baseq=${params.basequal} \\
        --qualityScore=${params.qualityscore} \\
        --maxNumberOfMismatchesInRead=${params.mismatch_threshold} \\
        --altBaseQualFile=${prefix}_alternative_allele_base_qualities.txt \\
        --refBaseQualFile=${prefix}_reference_allele_base_qualities.txt \\
        --altBasePositionsFile=${prefix}_alternative_allele_read_positions.txt \\
        --refBasePositionsFile=${prefix}_reference_allele_read_positions.txt \\
        --referenceFile=$fasta | \\
            confidenceAnnotation_SNVs.py $controlflag\\
                -i - \\
                $confoptions \\
                -a 0 \\
                --gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} \\
                --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} \\
                --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} \\
                --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} \\
                --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf} \\
                -f ${prefix}_somatic_snvs_for_bias.vcf > ${prefix}_peoverlap.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """     
}