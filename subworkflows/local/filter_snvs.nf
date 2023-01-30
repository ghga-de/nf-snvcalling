//
// FILTER SNVs: Filtering options
//

params.options = [:]

include { FILTER_BY_CRIT       } from '../../modules/local/filter_by_crit.nf'       addParams( options: params.options )
//include { INDEL_EXTRACTION     } from '../../modules/local/indel_extraction.nf'     addParams( options: params.options )
//include { VISUALIZE            } from '../../modules/local/visualize.nf'            addParams( options: params.options )
//include { INDEL_JSON           } from '../../modules/local/indel_json.nf'           addParams( options: params.options )

workflow FILTER_SNVS {
    take:
    vcf_ch        // channel: [val(meta), vcf]
    ref           // reference channel [ref.fa, ref.fa.fai]

    main:

    versions=Channel.empty()

    //
    // MODULE: FILTER_BY_CRIT
    //
    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
    vcf_ch
    )
    versions = versions.mix(FILTER_BY_CRIT.out.versions)

    emit:
    versions
}