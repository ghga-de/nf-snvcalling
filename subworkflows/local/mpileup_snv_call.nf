//
// SNVCALL: RUN BCFTOOLS MPILEUP by intervals
//

params.options = [:]

include { BCFTOOLS_MPILEUP   } from '../../modules/nf-core/modules/bcftools/mpileup/main'  addParams( options: params.options )
include { FILTER_STRAND_BIAS } from '../../modules/local/filter_strand_bias.nf'            addParams( options: params.options )
include { MPILEUP_COMPARE    } from '../../modules/local/mpileup_compare.nf'               addParams( options: params.options )
include { FILE_CONCATENATOR  } from '../../modules/local/file_concatenator.nf'             addParams( options: params.options )


workflow MPILEUP_SNV_CALL {
    take:
    sample_ch     // channel: [val(meta), tumor,tumor_bai, control, control_bai, tumorname, controlname]
    ref
    intervals

    main:
    versions = Channel.empty()

    // Combine intervals with samples to create 'interval x sample' number of parallel run
    sample_ch
        .combine(intervals)
        .set { combined_inputs }

    // RUN bcftools mpileup and bcftools call to call variants. This process is scattered by chr intervals
    BCFTOOLS_MPILEUP (
        combined_inputs, ref
    )
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    BCFTOOLS_MPILEUP.out.vcf
                    .join(BCFTOOLS_MPILEUP.out.stat)
                    .filter{meta, intervals, vcf, tbi, stat -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stat) > 0 }
                    .set{ch_vcf}
    ch_vcf.view()

    // RUN seqContext_annotator.pl and filterVcfForBias.py
    FILTER_STRAND_BIAS(
        BCFTOOLS_MPILEUP.out.vcf, ref
    )
    versions = versions.mix(FILTER_STRAND_BIAS.out.versions) 

    // RUN bcftools mpileup and vcf_pileup_compare_allin1_basecount.pl to compare germline variants.
    // This process only applies of there is control and runCompareGermline is true
    MPILEUP_COMPARE(
        FILTER_STRAND_BIAS.out.vcf, ref
    )
    versions = versions.mix(MPILEUP_COMPARE.out.versions) 

    // Group interval VCF files according to meta
    MPILEUP_COMPARE
        .out
        .vcf
        .groupTuple()
        .set { combined_vcf }

    // MERGE interval VCF files
    // RUN: headeredFileConcatenator.pl
    
    FILE_CONCATENATOR(
        combined_vcf
    )
    versions = versions.mix(FILE_CONCATENATOR.out.versions) 
    vcf_ch=FILE_CONCATENATOR.out.vcf

    emit:
    vcf_ch
    versions
}
