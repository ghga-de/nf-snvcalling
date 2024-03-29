process GET_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/odcf_mpileupsnvcalling:v0':'kubran/odcf_mpileupsnvcalling:v0' }"

    input:
    tuple val(meta), path(tumor), path(tumor_bai), path(control),  path(control_bai)

    output:
    path "contigs.bed"   , emit: contigs
    path "versions.yml"  , emit: versions

    script: 
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    
    if (params.runcontigs == 'ALT_HLA')
    {
        """
        touch ALT_contigs.bed
        touch HLA_contigs.bed

        samtools view -H $meta.tumor_bam | grep -P "_alt\\tLN" | sed -e 's/@SQ\\tSN://' -e 's/\\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 > ALT_contigs.bed
        samtools view -H $meta.tumor_bam | grep -P "SN:HLA" | sed -e 's/@SQ\\tSN://' -e 's/\\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 > HLA_contigs.bed

        cat ALT_contigs.bed HLA_contigs.bed > contigs.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS

        """
    }
    else {
        """
        samtools view -H $meta.tumor_bam | grep -P "SN:"  | sed -e 's/@SQ\\tSN://' -e 's/\tLN:/\\t0\\t/' -e 's/\\tAH.*//' | sort -V -k1,1 | cut -f 1 | tail -n +25 > contigs.bed
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools 2>&1) | sed -e 's/.*Version: //; s/ Usage.*//')
        END_VERSIONS
        """
    }
}
