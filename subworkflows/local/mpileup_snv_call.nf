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
    ref           // channel: [path(fasta), path(fai)]
    intervals     // channel: [[chr, region], [chr, region], ...]

    main:
    versions = Channel.empty()

    // Combine intervals with samples to create 'interval x sample' number of parallel run
    sample_ch
        .combine(intervals)
        .set { combined_inputs }

    //
    // MODULE:BCFTOOLS_MPILEUP 
    //
    // RUN bcftools mpileup and bcftools call to call variants. This process is scattered by chr intervals
    BCFTOOLS_MPILEUP (
        combined_inputs, ref
    )
    versions = versions.mix(BCFTOOLS_MPILEUP.out.versions)

    // filter VCFs if there is no variant
    BCFTOOLS_MPILEUP.out.vcf
                    .join(BCFTOOLS_MPILEUP.out.intervals)
                    .join(BCFTOOLS_MPILEUP.out.stats)
                    .filter{meta, vcf, intervals, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
                    .set{ch_vcf_stats}
    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, vcf]} 
        .set {ch_vcf}

    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, intervals]} 
        .set {ch_intervals} 

    //
    // MODULE:FILTER_STRAND_BIAS 
    //
    // RUN seqContext_annotator.pl and filterVcfForBias.py
    FILTER_STRAND_BIAS(
        ch_vcf.join(ch_intervals, by: [0]), ref
    )
    versions = versions.mix(FILTER_STRAND_BIAS.out.versions) 

    // filter VCFs if there is no variant after bias filtration
    FILTER_STRAND_BIAS.out.vcf
                    .join(FILTER_STRAND_BIAS.out.intervals)
                    .join(FILTER_STRAND_BIAS.out.stats)
                    .filter{meta, vcf, intervals, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
                    .set{ch_vcf_stats}
    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, vcf]} 
        .set {ch_vcf}

    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, intervals]} 
        .set {ch_intervals} 
    //
    // MODULE:MPILEUP_COMPARE 
    //
    // RUN bcftools mpileup and vcf_pileup_compare_allin1_basecount.pl to compare germline variants.
    // This process only applies of there is control and runCompareGermline is true
    MPILEUP_COMPARE(
        ch_vcf.join(ch_intervals, by: [0]), ref
    )
    versions = versions.mix(MPILEUP_COMPARE.out.versions) 

    // Group interval VCF files according to meta
    MPILEUP_COMPARE
        .out
        .vcf
        .groupTuple()
        .set { combined_vcf }

    //
    // MODULE: FILE_CONCATENATOR 
    //
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
