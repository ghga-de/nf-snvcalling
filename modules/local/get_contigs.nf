process GET_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)
    tuple val(meta2), path(contig_file)

    output:
    tuple val(meta), path("contigs.bed") , emit: contigs
    path "versions.yml"                  , emit: versions

    script: 
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    if (params.contig_file)
    {
        """
        cp $contig_file contigs.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """
    }
    else{
        if (params.runcontigs == 'ALT_HLA')
        {
            """
            if [ -n "\$(samtools view -H $tumor | grep -P "_alt\\tLN")" ]; then
                samtools view -H $tumor | grep -P "_alt\\tLN" | sed -e 's/@SQ\\tSN://' -e 's/\\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 > ALT_contigs.bed
            else
                touch ALT_contigs.bed
            fi
            if [ -n "\$(samtools view -H $tumor | grep -P "SN:HLA")" ]; then
                samtools view -H $tumor | grep -P "SN:HLA" | sed -e 's/@SQ\\tSN://' -e 's/\\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 > HLA_contigs.bed
            else
                touch HLA_contigs.bed
            fi
            cat ALT_contigs.bed HLA_contigs.bed > contigs.bed

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
            END_VERSIONS

            """
        }
        else {
            """
            touch contigs.bed
            if [ -n "\$(samtools view -H $tumor | grep -P "SN:") "]; then
                samtools view -H $tumor | grep -P "SN:"  | sed -e 's/@SQ\\tSN://' -e 's/\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 | tail -n +25 > contigs.bed
            else
                touch contigs.bed
            fi
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
            END_VERSIONS
            """
        }

    }
}