// Applied only for no-control cases
process FILTER_BY_CRIT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(vcfgz), path(vcf_tbi)

    output:
    tuple val(meta), path('*Filtered.vcf.gz'),  path('*Filtered.vcf.gz.tbi')   , emit: vcf
    path  "versions.yml"                                                       , emit: versions
    path  "*_postFilter_criteria.txt"                                          , optional: true 

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    // recurrance is not implemented!!
    def filter_values = [( params.gnomad_genomes && params.crit_gnomad_genomes_maxmaf ) ? "GNOMAD_GENOMES AF $params.crit_gnomad_genomes_maxmaf+": "",
                        ( params.gnomad_exomes && params.crit_gnomad_exomes_maxmaf) ? "GNOMAD_EXOMES AF $params.crit_gnomad_exomes_maxmaf+": "",
                        ( params.k_genome && params.crit_1kgenomes_maxmaf) ? "1K_GENOMES EUR_AF $params.crit_1kgenomes_maxmaf+": "",
                        ( params.dbsnp_snv && params.filter_non_clinic)? "DBSNP CLN,COMMON nonexist,exist": "", 
                        ( params.local_control_wgs && params.local_control_wes && params.crit_localcontrol_maxmaf) ? "LocalControlAF_WGS AF $params.crit_localcontrol_maxmaf+ LocalControlAF_WES AF $params.crit_localcontrol_maxmaf+": "",
                        ].join(' ').trim() 

// Filter variants only if there is no control, else do noting
    if (meta.iscontrol == '1') {
        """
        mv $vcfgz ${prefix}_noFiltered.vcf.gz
        tabix ${prefix}_noFiltered.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
        END_VERSIONS
        """
    }
    else
        {
        if (filter_values) 
            {
            """
            vcf_filter_by_crit.py $vcfgz ${prefix}_postFiltered.vcf $filter_values
            bgzip -c ${prefix}_postFiltered.vcf > ${prefix}_postFiltered.vcf.gz
            tabix ${prefix}_postFiltered.vcf.gz

            echo "$filter_values" > ${prefix}_postFilter_criteria.txt

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python2 --version 2>&1 | sed 's/Python //g')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            END_VERSIONS
            """
            }
        else
        {
            """
            mv $vcfgz ${prefix}_noFiltered.vcf.gz
            tabix ${prefix}_noFiltered.vcf.gz
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            END_VERSIONS
            """ 
        }

    }

}