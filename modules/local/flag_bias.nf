
process FLAG_BIAS {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), path(vcf), path(sequence_error_matrix), path(sequencing_error_matrix)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*.flagbias.vcf'),    emit: vcf
    tuple val(meta), path('*.flagbiastmp.vcf'), emit: vcftmp
    path  "versions.yml"                   ,    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag = meta.iscontrol == "1" ? "true" : "false"
    def ref_hg38   = params.ref_type == "hg38" ? "--refgenome GRCh38 ftp://ftp.sanger.ac.uk/pub/cancer/dockstore/human/GRCh38_hla_decoy_ebv/core_ref_GRCh38_hla_decoy_ebv.tar.gz": "" 
    def ref_spec   = params.ref_type == "hg38" ? "--gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf}" : ""

    """
    filterVcfForBias.py \\
        --vcfFile=$vcf \\
        --referenceFile=$fasta
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
        --bias_matrixSeqFile=${filenameBiasMatrixSeqFile} \\
        --bias_matrixSeqingFile=${filenameBiasMatrixSeqingFile}
        --vcfFileFlagged=${prefix}.flagged.vcf | \\
            confidenceAnnotation_SNVs.py $controlflag -i - ${params.confidenceoptions} $ref_hg38 -a 1 \\
                -f snv_${prefix}.flagbiastmp.vcf $ref_spec > snv_${prefix}.flagbias.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}