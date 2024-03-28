process CONVERT_TO_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(input), path(raw_vcf)
    path(config)

    output:
    tuple val(meta), path("*.std.vcf") ,  emit: std_vcf
    path  "versions.yml"               ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf_name = input.getExtension() == "gz" ? input.getBaseName() : input.getName()
    vcf_name = vcf_name.take(vcf_name.size() - 3)
    def header = raw_vcf ? "-r $raw_vcf" : "-r False"

    """
    convertToStdVCF.py -i $input \\
        $header \\
        -s $prefix \\
        -w False \\
        -c $config \\
        -o ${vcf_name}std.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python2 --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
