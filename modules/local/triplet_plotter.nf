// This process is not running since BaseQualitiesForChosenPositions.1000000.1.txt.gz file which needs to be in aligment fle is missing!
// This process will be asked!

process TRIPLET_PLOTTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v8':'kubran/odcf_snvcalling:v8' }"

    input:
    tuple val(meta), path(somaticvcf), path(vcf), path(index),  path(altbasequal), path(refbasequal), path(altreadpos), path(refreadpos)
    val(plotBackgroundBaseScoreDistribution)
    val(title)

    output:
    path "*withMAF.vcf"             , emit: withmaf_vcf
    path "*filteredAltMedian*.vcf"  , emit: filtered_vcf
    path "*.pdf"                    , emit: pdf
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def force = params.rerunfiltering ? '1': '0'
    def median_filter_threshold = params.rerunfiltering  ? "${params.median_filter_threshold}" : '-1'
    def skipplots    = params.rerunfiltering  ? '1': '0'
    def rerun_suffix = params.rerunfiltering  ? '1': '0'

    """
    tripletbased_BQdistribution_runner.sh \\
        -p ${prefix} \\
        -i $somaticvcf \\
        -t $meta.tumor_bam \\
        -v ${params.basequal} \\
        -b $plotBackgroundBaseScoreDistribution \\
        -o ${prefix}_tripletSpecific_base_score_distribution${rerun_suffix} \\
        -r $force \\
        -m $median_filter_threshold \\
        -w "${prefix}: ${title}" \\
        -sp $skipplots  \\
        -rb $refbasequal \\
        -ab $altbasequal \\
        -ar $altreadpos \\
        -rr $refreadpos

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
        python: \$(python2.7 --version | sed 's/Python //g')
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools //; s/Using.*\$//') 
    END_VERSIONS
    """
}
