// Applied only for no-control cases
process FILTER_BY_CRIT {
    tag "$meta.id"
    label 'process_medium'

    conda     (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'docker://kubran/odcf_snvcalling:v7':'kubran/odcf_snvcalling:v7' }"

    input:
    tuple val(meta), file(vcfgz), file(vcf_tbi)

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
    def filter_values = [ (params.filter_exac && params.crit_exac_maxmaf) ? "ExAC AF $params.crit_exac_maxmaf+": "",
                        (params.filter_evs && params.crit_evs_maxmaf) ? "EVS MAF $params.crit_evs_maxmaf+": "",
                        (params.filter_gnomad_genomes && params.crit_gnomad_genomes_maxmaf) ? "GNOMAD_GENOMES AF $params.crit_gnomad_genomes_maxmaf+": "",
                        (params.filter_gnomad_exomes && params.crit_gnomad_exomes_maxmaf) ? "GNOMAD_EXOMES AF $params.crit_gnomad_exomes_maxmaf+": "",
                        (params.filter_1kgenomes && params.crit_1kgenomes_maxmaf) ? "1K_GENOMES AF $params.crit_1kgenomes_maxmaf+": "",
                        params.filter_non_clinic ? "DBSNP CLN,COMMON nonexist,exist": "", 
                        (params.filter_localcontrol && params.crit_localcontrol_maxmaf) ? "LocalControlAF_WGS AF $params.crit_localcontrol_maxmaf+": "",
                        (params.filter_localcontrol && params.crit_localcontrol_maxmaf) ? "LocalControlAF_WES AF $params.crit_localcontrol_maxmaf+": ""              
                        ].join(' ').trim() 

// Filter variants only if there is no control, else do noting
    if (meta.iscontrol == '1') {
        """
        mv $vcfgz ${prefix}_noFiltered.vcf.gz
        tabix ${prefix}_noFiltered.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
            gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
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
                python: \$(python --version | sed 's/Python //g')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
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
                python: \$(python --version | sed 's/Python //g')
                tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
                gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
            END_VERSIONS
            """ 
        }

    }

}