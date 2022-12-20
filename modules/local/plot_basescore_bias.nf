process PLOT_BASESCORE_BIAS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"

    input:
    tuple val(meta), path(sequencing_error_matrix), path(sequence_error_matrix)

    output:
    path  "versions.yml"                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (params.runArtifactFilter){
        """
        plotBaseScoreDistibution.R \\
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
            --ratioVcf=${paramsrVcf} \\
            --bias_matrixSeqFile=${filenameBiasMatrixSeqFile} 
            --bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile} --vcfFileFlagged="/dev/stdout" | \
	${PYPY_OR_PYTHON_BINARY} -u ${TOOL_CONFIDENCE_ANNOTATION} ${noControlFlag} -i - ${CONFIDENCE_OPTS} -a 1 -f ${filenameSomaticSNVsTmp} \
        --gnomAD_WGS_maxMAF=${CRIT_GNOMAD_GENOMES_maxMAF} --gnomAD_WES_maxMAF=${CRIT_GNOMAD_EXOMES_maxMAF} --localControl_WGS_maxMAF=${CRIT_LOCALCONTROL_maxMAF} --localControl_WES_maxMAF=${CRIT_LOCALCONTROL_maxMAF} --1000genome_maxMAF=${CRIT_1KGENOMES_maxMAF} > ${filenameSNVVCFTemp}.tmp

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            perl: v5.28.1
            python: \$(python2.7 --version | sed 's/Python //g')
            bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
        END_VERSIONS
        """
    }
}
