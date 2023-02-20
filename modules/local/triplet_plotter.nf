// This process will be asked!

process TRIPLET_PLOTTER {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_snvcalling:v10':'kubran/odcf_snvcalling:v10' }"

    input:
    tuple val(meta), path(somaticvcf), path(altbasequal), path(refbasequal), path(altreadpos), path(refreadpos)
    val(title)

    output:
    tuple val(meta), path("*.withMAF.vcf")             , emit: withmaf_vcf
    tuple val(meta), path("*filteredAltMedian*.vcf")  , emit: filtered_vcf
    tuple val(meta), path("*.pdf")                    , emit: plot
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def force = params.rerunfiltering ? '1': '0'
    def median_filter_threshold = params.rerunfiltering  ? "${params.median_filter_threshold}" : '-1'
    def skipplots    = params.rerunfiltering  ? '1': '0'
    def rerun_suffix = params.rerunfiltering  ? '1': "''"

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
        -R $force \\
        -c 1 \\
        -f $median_filter_threshold \\
        -s \$SEQUENCE_CONTEXT_COLUMN_INDEX \\
        --MAFColumnIndex \$MAF_COLUMN_INDEX \\
        -i 1 \\
        -b "0" \\
        --skipPlots $skipplots \\
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
