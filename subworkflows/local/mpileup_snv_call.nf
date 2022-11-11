//
// SNVCALL: RUN P
//

params.options = [:]

include { BCFTOOLS_MPILEUP   } from '../../modules/local/bcftools_mpileup.nf'         addParams( options: params.options )
include { FILTER_STRAND_BIAS } from '../../modules/local/filter_strand_bias.nf'       addParams( options: params.options )
include { MPILEUP_COMPARE    } from '../../modules/local/mpileup_compare.nf'          addParams( options: params.options )

workflow MPILEUP_SNV_CALL {
    take:
    sample_ch     // channel: [val(meta), tumor,tumor_bai, control, control_bai, tumorname, controlname]
    ref
    intervals

    main:
    versions = Channel.empty()
    sample_ch
        .combine(intervals)
        .set { combined_inputs }

    // RUN bcftools mpileup and bcftools call to call variants. This process is scattered by chr intervals
    BCFTOOLS_MPILEUP (
        combined_inputs, ref
    )
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    // RUN seqContext_annotator.pl and filterVcfForBias.py
    //FILTER_STRAND_BIAS(
    //    BCFTOOLS_MPILEUP.out.vcf, ref
    //)
    //versions = versions.mix(FILTER_STRAND_BIAS.out.versions) 

    // RUN bcftools mpileup and vcf_pileup_compare_allin1_basecount.pl to compare germline variants.
    // This process only applies of there is control and runCompareGermline is true

    //MPILEUP_COMPARE(
    //    FILTER_STRAND_BIAS.out.vcf, 
    //)

    // 


    // AT THE END CHR results will be merged
    // RUN: headeredFileConcatenator.pl
    

    emit:
    version
}
