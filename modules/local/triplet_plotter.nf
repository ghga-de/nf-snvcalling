
process TRIPLET_PLOTTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"

    input:
    tuple val(meta), path(somaticvcf), path(altbasequal), path(refbasequal), path(altreadpos), path(refreadpos)
    val(title)
    val(rerun)

    output:
    tuple val(meta), path("*withMAF_filtered*.vcf")   , emit: vcf    , optional: true
    tuple val(meta), path("*_combined.pdf")           , emit: plot_1 , optional: true
    tuple val(meta), path("*_CoV.pdf")                , emit: plot_2 , optional: true
    tuple val(meta), path("*_CHROMcolored.pdf")       , emit: plot_3 , optional: true
    tuple val(meta), path("*_Q50.pdf")                , emit: plot_4 , optional: true
    tuple val(meta), path("*_Q60.pdf")                , emit: plot_5 , optional: true
    tuple val(meta), path("*_Q70.pdf")                , emit: plot_6 , optional: true
    tuple val(meta), path("*_Q80.pdf")                , emit: plot_7 , optional: true
    tuple val(meta), path("*_VAFcolored.pdf")         , emit: plot_8 , optional: true
    tuple val(meta), path("*withMAF.vcf")             , emit: withmaf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def mft          = rerun == 1  ? "${params.median_filter_threshold}" : '-1'
    def skip_force   = rerun == 1  ? '1': '0'
    def rerun_suffix = rerun ? "_filteredAltMedian${params.median_filter_threshold}": ""     

    """
    cat $somaticvcf | perl -ne 'chomp; my \$line=\$_; if (/DP4=(\\d+),(\\d+),(\\d+),(\\d+);/) {my \$fR=\$1; my \$rR=\$2; my \$fA=\$3; my \$rA=\$4; my \$MAF=(\$fA+\$rA)/(\$fR+\$rR+\$fA+\$rA); print "\$line\\t\$MAF\\n";} else { if (/^#CHROM/) { print "\$line\\tMAF\\n";} else {print "\$line\\n";} };' >${prefix}.withMAF.vcf

    SEQUENCE_CONTEXT_COLUMN_INDEX=`cat ${prefix}.withMAF.vcf | grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, \$_); my \$columnIndex = first { \$colnames[\$_] eq "SEQUENCE_CONTEXT"} 0..\$#colnames; \$columnIndex += 1; print "\$columnIndex\\n";'`
    echo \$SEQUENCE_CONTEXT_COLUMN_INDEX
    MAF_COLUMN_INDEX=`cat ${prefix}.withMAF.vcf| grep -v '^##' | grep '^#' | perl -ne 'use List::Util qw(first); chomp; my @colnames = split(/\t/, \$_); my \$columnIndex = first { \$colnames[\$_] eq "MAF"} 0..\$#colnames; \$columnIndex += 1; print "\$columnIndex\\n";'`
    echo \$MAF_COLUMN_INDEX

    tripletBased_BQDistribution_plotter.R -v ${prefix}.withMAF.vcf \\
        -a "" \\
        -p $prefix \\
        -t "${prefix} ${title}" \\
        -o ${prefix}_tripletSpecific_base_score_distribution${rerun_suffix} \\
        -R $skip_force \\
        -c 1 \\
        -f $mft \\
        -s \$SEQUENCE_CONTEXT_COLUMN_INDEX \\
        --MAFColumnIndex \$MAF_COLUMN_INDEX \\
        -i 1 \\
        -b 0 \\
        --skipPlots $skip_force \\
        --refBaseQual $refbasequal \\
        --altBaseQual $altbasequal \\
        --altReadPos $altreadpos \\
        --refReadPos $refreadpos

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        perl: \$(echo \$(perl --version 2>&1) | sed 's/.*v\\(.*\\)) built.*/\\1/')
    END_VERSIONS
    """
}
