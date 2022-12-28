//# create BaseScore FIFOs and their consumer processes (zip and write to target file)
//# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}

process POST_PROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*.postprocessed.vcf.gz'), path('*.postprocessed.vcf.gz.tbi') , emit: vcf
    path  "versions.yml"                                                                 , emit: versions
    path ('*.txt')     
    path ('*.pdf')
    path ('*.vcf')

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag = meta.iscontrol == "1" ? "true" : "false"
    def ref_spec   = params.ref_type == "hg38" ? "'--gnomAD_WGS_maxMAF=${params.crit_gnomad_genomes_maxmaf} --gnomAD_WES_maxMAF=${params.crit_gnomad_exomes_maxmaf} --localControl_WGS_maxMAF=${params.crit_localcontrol_maxmaf} --localControl_WES_maxMAF=${params.crit_localcontrol_maxmaf} --1000genome_maxMAF=${params.crit_1kgenomes_maxmaf}'" : ""
    """
    post_process.sh \\
        -i $vcfgz \\
        -p $prefix \\
        -t $params.runArtifactFilter \\
        -q $params.generateExtendedQcPlots \\
        -c $controlflag \\
        -tb $meta.tumor_bam \\
        -ofpe '--mapq=${params.mapequal} --baseq=${params.basequal} --qualityScore=${params.qualityScore} --maxNumberOfMismatchesInRead=${params.mismatchesinread}' \\
        -ovb '--numReads=${params.nReads} --numMuts=${params.nMuts} --biasPValThreshold=${params.biasPValThreshold} --biasRatioThreshold=${params.biasRatioThreshold} --biasRatioMinimum=${params.biasRatioMinimum} --maxNumOppositeReadsSequencingWeakBias=${params.maxNumOppositeReadsSequencingWeakBias} --maxNumOppositeReadsSequenceWeakBias=${params.maxNumOppositeReadsSequenceWeakBias} --maxNumOppositeReadsSequencingStrongBias=${params.maxNumOppositeReadsSequencingStrongBias} --maxNumOppositeReadsSequenceStrongBias=${params.maxNumOppositeReadsSequenceStrongBias} --ratioVcf=${params.rVcf}' \\
        -r $fasta \\
        -oc '$args' \\
        -rs '$ref_spec' \\
        -bq $params.basequal

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}