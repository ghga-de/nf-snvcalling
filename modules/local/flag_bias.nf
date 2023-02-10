process FLAG_BIAS {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(vcf), file(sequence_error_matrix), file(sequencing_error_matrix)
    tuple path(fasta), path(fai)
    val(round)
    val(bias_round)

    output:
    tuple val(meta), path("*.${bias_round}.flagged.vcf"),       emit: flagged
    tuple val(meta), path("*.${bias_round}.flagbias.vcf"),      emit: vcf
    tuple val(meta), path("*.${bias_round}.flagbiastmp.vcf"),   emit: vcftmp
    tuple val(meta), path('*.txt'),                             emit: bias_matrixs
    path  "versions.yml" ,                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag = meta.iscontrol == "1" ? "" : "--nocontrol"
    def confoptions   = params.ref_type == "hg38" ? "${params.confidenceoptions} --refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "${params.confidenceoptions}" 
    
    """
    filterVcfForBias.py \\
        --vcfFile=$vcf \\
        --referenceFile=$fasta \\
        --sequence_specificFile=$sequence_error_matrix \\
        --sequencing_specificFile=$sequencing_error_matrix \\
        --numReads=${params.nReads} \\
        --numMuts=${params.nMuts} \\
        --biasPValThreshold=${params.biasPValThreshold} \\
        --biasRatioThreshold=${params.biasRatioThreshold} \\
        --biasRatioMinimum=${params.biasRatioMinimum} \\
        --maxNumOppositeReadsSequencingWeakBias=${params.maxNumOppositeReadsSequencingWeakBias} \\
        --maxNumOppositeReadsSequenceWeakBias=${params.maxNumOppositeReadsSequenceWeakBias} \\
        --maxNumOppositeReadsSequencingStrongBias=${params.maxNumOppositeReadsSequencingStrongBias} \\
        --maxNumOppositeReadsSequenceStrongBias=${params.maxNumOppositeReadsSequenceStrongBias} \\
        --ratioVcf=${params.rVcf} \\
        --bias_matrixSeqFile=${prefix}_sequence_specific_bias_matrix_${bias_round}.txt \\
        --bias_matrixSeqingFile=${prefix}_sequencing_specific_bias_matrix_${bias_round}.txt \\
        --vcfFileFlagged=${prefix}.${bias_round}.flagged.vcf

        cat < ${prefix}.${bias_round}.flagged.vcf | confidenceAnnotation_SNVs.py \\
            $controlflag \\
            -i - \\
            $confoptions \\
            -a $round \\
            --gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} \\
            --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} \\
            --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} \\
            --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf} \\
            -f snv_${prefix}.${bias_round}.flagbiastmp.vcf > snv_${prefix}.${bias_round}.flagbias.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2.7 --version | sed 's/Python //g')
    END_VERSIONS
    """
}