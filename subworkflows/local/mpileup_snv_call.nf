//
// SNVCALL: RUN BCFTOOLS MPILEUP by intervals
//

params.options = [:]

include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_1 } from '../../modules/nf-core/modules/bcftools/mpileup/main'  addParams( options: params.options )
include { BCFTOOLS_MPILEUP as BCFTOOLS_MPILEUP_2 } from '../../modules/nf-core/modules/bcftools/mpileup/main'  addParams( options: params.options )
include { SEQ_CONTEXT_ANNOTATOR } from '../../modules/local/seq_context_annotator.nf'         addParams( options: params.options )
include { MPILEUP_COMPARE       } from '../../modules/local/mpileup_compare.nf'               addParams( options: params.options )
include { FILE_CONCATENATOR     } from '../../modules/local/file_concatenator.nf'             addParams( options: params.options )


workflow MPILEUP_SNV_CALL {
    take:
    sample_ch     // channel: [val(meta), tumor,tumor_bai, control, control_bai, tumorname, controlname]
    ref           // channel: [path(fasta), path(fai)]
    intervals     // channel: [[chr, region], [chr, region], ...]
    contigs       // channel: [path(hla)]

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
    combined_inputs = combined_inputs.map {it -> tuple( it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7], [])} 
    BCFTOOLS_MPILEUP_1 (
        combined_inputs, ref
    )
    versions = versions.mix(BCFTOOLS_MPILEUP_1.out.versions)

    // filter VCFs if there is no variant
    BCFTOOLS_MPILEUP_1.out.vcf
                    .join(BCFTOOLS_MPILEUP_1.out.intervals)
                    .join(BCFTOOLS_MPILEUP_1.out.stats)
                    .filter{meta, vcf, intervals, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
                    .set{ch_vcf_stats}
    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, vcf]} 
        .set {ch_vcf_1}
    ch_vcf_stats
        .map { meta, vcf, intervals, stats -> [meta, intervals]} 
        .set {ch_intervals_1}
    // Collect VCF files and intervals
    ch_vcf = Channel.empty()
    ch_intervals = Channel.empty()
    ch_vcf = ch_vcf.mix(ch_vcf_1)
    ch_intervals = ch_intervals.mix(ch_intervals_1)

    // If reference is hg38, run HLA and ALT contigs as seperate runs
    if (params.ref_type == 'hg38'){

        // Generate input set for HLA and ALT runs
        hla_alt_inputs = sample_ch.combine(contigs)
        hla_alt_inputs = hla_alt_inputs.map {it -> tuple( it[0], it[1], it[2], it[3], it[4], it[5], it[6], [], it[7])} 
        BCFTOOLS_MPILEUP_2 (
            hla_alt_inputs, ref
        )
        BCFTOOLS_MPILEUP_2.out.vcf
                    .join(BCFTOOLS_MPILEUP_2.out.intervals)
                    .join(BCFTOOLS_MPILEUP_2.out.stats)
                    .filter{meta, vcf, intervals, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
                    .set{ch_vcf_stats_hla_alt}
        ch_vcf_stats_hla_alt
            .map { meta, vcf, intervals, stats -> [meta, vcf]} 
            .set {ch_vcf_2}
        ch_vcf = ch_vcf.mix(ch_vcf_2)

        ch_vcf_stats_hla_alt
            .map { meta, vcf, intervals, stats -> [meta, intervals]} 
            .set {ch_intervals_2}
        ch_intervals = ch_intervals.mix(ch_intervals_2)
    }

    //
    // MODULE:SEQ_CONTEXT_ANNOTATOR
    //
    // RUN seqContext_annotator.pl and filterVcfForBias.py
    SEQ_CONTEXT_ANNOTATOR(
        ch_vcf.join(ch_intervals, by: [0]), ref
    )
    //versions = versions.mix(SEQ_CONTEXT_ANNOTATOR.out.versions) 

    // filter VCFs if there is no variant after bias filtration
    SEQ_CONTEXT_ANNOTATOR.out.vcf
                    .join(SEQ_CONTEXT_ANNOTATOR.out.intervals)
                    .join(SEQ_CONTEXT_ANNOTATOR.out.stats)
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
    //versions = versions.mix(MPILEUP_COMPARE.out.versions) 

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
    //versions = versions.mix(FILE_CONCATENATOR.out.versions) 
    vcf_ch=FILE_CONCATENATOR.out.vcf

    emit:
    vcf_ch
    versions
}
