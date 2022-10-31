//
// SNVCALL: RUN P
//

params.options = [:]

include { BCFTOOLS_MPILEUP   } from '../../modules/nf-core/modules/bcftools/mpileup/main.nf'         addParams( options: params.options )
//include { CHECK_IF_CORRUPTED} from '../../modules/local/check_if_corrupted.nf'     addParams( options: params.options )

workflow SNV_CALL {
    take:
    sample_ch     // channel: [val(meta), tumor,tumor_bai, control, control_bai, tumorname, controlname]
    ref
    intervals

    main:
    BCFTOOLS_MPILEUP (
    sample_ch, ref, intervals
    )
    version = BCFTOOLS_MPILEUP.out.versions

    emit:
    version
}
