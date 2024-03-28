//
// OUTPUT_STANDARD_VCF: OUTPUT_STANDARD_VCF
//

params.options = [:]
 
include { CONVERT_TO_VCF as CONVERT_TO_VCF_1 } from '../../modules/local/convert_to_vcf.nf'  addParams( options: params.options )   
include { CONVERT_TO_VCF as CONVERT_TO_VCF_2 } from '../../modules/local/convert_to_vcf.nf'  addParams( options: params.options )   
include { BCFTOOLS_SORT         } from '../../modules/nf-core/modules/bcftools/sort/main'    addParams( options: params.options )          
include { TABIX_BGZIPTABIX      } from '../../modules/nf-core/modules/tabix/bgziptabix/main' addParams( options: params.options )             

workflow OUTPUT_STANDARD_VCF {
    take:
    vcf_ch // channel: [val(meta), vcf, header]
    config
    
    main:

    versions=Channel.empty()
    //
    // MODULE: CONVERT_TO_VCF
    //
    CONVERT_TO_VCF_1(
        vcf_ch,
        config
    )
    versions = versions.mix(CONVERT_TO_VCF_1.out.versions)

    CONVERT_TO_VCF_2(
        CONVERT_TO_VCF_1.out.std_vcf.map{ it -> tuple( it[0], it[1], [] )},
        config
    )

    TABIX_BGZIPTABIX(
        CONVERT_TO_VCF_2.out.std_vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX.out.versions)

    BCFTOOLS_SORT(
        TABIX_BGZIPTABIX.out.gz_tbi
    )
    versions = versions.mix(BCFTOOLS_SORT.out.versions)

    emit:
    versions
}
