//# create BaseScore FIFOs and their consumer processes (zip and write to target file)
//# BaseScore FIFOS will be filled by ${TOOL_FILTER_PE_OVERLAP}

process FILTER_PEOVERLAP {
    tag "$meta.id"
    label 'process_high'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v2':'kubran/odcf_snvcalling:v2' }"
    
    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path('*.allelebasescore.vcf')  , emit: allelebasescore_vcf
    path  "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def controlflag = meta.iscontrol == "1" ? "true" : "false"

    """
    filter_peoverlap.sh -i $vcfgz \\
        -p $prefix \\
        -t $params.mismatchesinread \\
        -gnmax $params.crit_gnomad_genomes_maxmaf \\
        -gemax $params.crit_gnomad_exomes_maxmaf \\
        -lcmax $params.crit_localcontrol_maxmaf \\
        -kgmax $params.crit_1kgenomes_maxmaf \\
        -c $controlflag \\
        -tb $meta.tumor_bam \\
        -mq $params.mapqual \\
        -bq $params.basequal \\
        -r $fasta \\
        -o ${prefix}.allelebasescore.vcf
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS
    """
}