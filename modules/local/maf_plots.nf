process MAF_PLOTS {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"
    
    input:
    tuple val(meta), file(somatic_snv), file(maf_values), file(indbsnp)

    output:
    tuple val(meta), path('*.pdf')       , emit: plot
    tuple val(meta), path('*.json')      , emit: json  
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    """
    MAF_plots.r \\
        $maf_values \\
        `grep -v "^#" ${somatic_snv} | wc -l` \\
        ${prefix}_MAF_conf_${params.min_confidence_score}_to_10.pdf \\
        $prefix \\
        ` awk '{FS="\t"}{if(NR==2)print \$5}'	${indbsnp}`

    snvnum=`grep -v "^#" ${somatic_snv} | wc -l`
    snvindbSNP=` awk '{FS="\t"}{if(NR==2)print \$5}'	${indbsnp}`
    SNV_IN_DBSNP_RATIO=`echo -e "\$snvindbSNP\t\$snvnum" | perl -F -ne 'print \$F[0]/\$F[1];'`

    echo -e "{" >${prefix}_QC_values.json
    echo -e "\t\"snvnum\": ${snvnum:-NA}," >>${prefix}_QC_values.json
    echo -e "\t\"snvindbSNP\": ${snvindbSNP:-NA}," >>${prefix}_QC_values.json
    echo -e "\t\"snvInDbsnpRatio\": ${SNV_IN_DBSNP_RATIO:-NA}," >>${prefix}_QC_values.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: v5.28.1
    END_VERSIONS
    """
}