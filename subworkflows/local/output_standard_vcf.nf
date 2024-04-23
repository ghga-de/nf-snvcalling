//
// OUTPUT_STANDARD_VCF: OUTPUT_STANDARD_VCF
//

params.options = [:]
 
include { CONVERT_TO_VCF    } from '../../modules/local/convert_to_vcf.nf'               addParams( options: params.options )   
include { BCFTOOLS_SORT     } from '../../modules/nf-core/modules/bcftools/sort/main'    addParams( options: params.options )          
include { TABIX_BGZIPTABIX  } from '../../modules/nf-core/modules/tabix/bgziptabix/main' addParams( options: params.options )             

workflow OUTPUT_STANDARD_VCF {
    take:
    vcf_ch // channel: [val(meta), vcf, header]
    config
    
    main:

    versions=Channel.empty()
    //
    // MODULE: CONVERT_TO_VCF
    //
    // use convertToStdVCF.json to convert DKFZ files to standard VCF4.2
    CONVERT_TO_VCF(
        vcf_ch,
        config
    )
    versions = versions.mix(CONVERT_TO_VCF.out.versions)

    //
    // MODULE: TABIX_BGZIPTABIX
    //
    TABIX_BGZIPTABIX(
        CONVERT_TO_VCF.out.std_vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX.out.versions)

    //
    // MODULE: BCFTOOLS_SORT
    //
    BCFTOOLS_SORT(
        TABIX_BGZIPTABIX.out.gz_tbi
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    versions
}
