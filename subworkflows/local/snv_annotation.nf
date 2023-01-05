//
// SNV ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF           } from '../../modules/local/annotate_vcf.nf'            addParams( options: params.options )
include { ANNOVAR                } from '../../modules/local/annovar.nf'                 addParams( options: params.options )
include { SNV_RELIABILITY_PIPE   } from '../../modules/local/snv_reliability_pipe.nf'    addParams( options: params.options )
include { CONFIDENCE_ANNOTATION  } from '../../modules/local/confidence_annotation.nf'   addParams( options: params.options )
include { ANNOTATION_PIPES       } from '../../modules/local/annotation_pipes.nf'        addParams( options: params.options )
include { POST_PROCESS           } from '../../modules/local/post_process.nf'            addParams( options: params.options )
include { FILTER_PEOVERLAP       } from '../../modules/local/filter_peoverlap.nf'        addParams( options: params.options )


workflow SNV_ANNOTATION {
    take:
    vcf_ch               // channel: [val(meta), , vcf.gz, vcf.gz.tbi ,val(tumorname), val(controlname) ]
    ref                  // channel: [path(fasta), path(fai)]
    kgenome              // channel: [file.vcf.gz, file.vcf.gz.tbi]
    dbsnpsnv             // channel: [file.vcf.gz, file.vcf.gz.tbi]
    localcontrolwgs      // channel: [file.vcf.gz, file.vcf.gz.tbi]
    localcontrolwes      // channel: [file.vcf.gz, file.vcf.gz.tbi]
    gnomadgenomes        // channel: [file.vcf.gz, file.vcf.gz.tbi]
    gnomadexomes         // channel: [file.vcf.gz, file.vcf.gz.tbi]
    annodb               // channel: [table_annovar_dir]
    repeatmasker         // channel: [file.bed.gz, file.bed.gz.tbi]
    dacblacklist         // channel: [file.bed.gz, file.bed.gz.tbi]
    dukeexcluded         // channel: [file.bed.gz, file.bed.gz.tbi]
    hiseqdepth           // channel: [file.bed.gz, file.bed.gz.tbi]
    selfchain            // channel: [file.bed.gz, file.bed.gz.tbi]
    mapability           // channel: [file.bed.gz, file.bed.gz.tbi]
    simpletandemrepeats  // channel: [file.bed.gz, file.bed.gz.tbi]
    enchangers           // channel: [file.bed.gz, file.bed.gz.tbi]
    cpgislands           // channel: [file.bed.gz, file.bed.gz.tbi]
    tfbscons             // channel: [file.bed.gz, file.bed.gz.tbi]
    encode_dnase         // channel: [file.bed.gz, file.bed.gz.tbi]
    mirnas_snornas       // channel: [file.bed.gz, file.bed.gz.tbi]
    cosmic               // channel: [file.bed.gz, file.bed.gz.tbi]
    mirbase              // channel: [file.bed.gz, file.bed.gz.tbi]
    mir_targets          // channel: [file.bed.gz, file.bed.gz.tbi]
    cgi_mountains        // channel: [file.bed.gz, file.bed.gz.tbi]
    phastconselem        // channel: [file.bed.gz, file.bed.gz.tbi]
    encode_tfbs          // channel: [file.bed.gz, file.bed.gz.tbi]
    mirnas_sncrnas       // channel: [file.bed.gz, file.bed.gz.tbi] 
    chr_prefix           // val channel: [prefix]

    main:

    versions=Channel.empty()
    logs=Channel.empty() 

    //
    // MODULE: ANNOTATE_VCF
    //
    // RUN annotate_vcf.pl: Uses various databases (all mandatory exept recurrance) to annotate variants
    ANNOTATE_VCF (
        vcf_ch, kgenome, dbsnpsnv, localcontrolwgs,
        localcontrolwes, gnomadgenomes, gnomadexomes, chr_prefix
    )
    versions  = versions.mix(ANNOTATE_VCF.out.versions)

    // RUN annovar, processAnnovarOutput.pl and newCols2vcf.pl: annovar annotates and classifies the variants, 
    // perl scripts re-creates vcfs.

    ch_vcf = ANNOTATE_VCF.out.unziped_vcf
    input_ch = ch_vcf.join(ANNOTATE_VCF.out.forannovar)

    // Only if RunGeneAnnovar is true, run annovar
    //
    // MODULE: ANNOVAR
    //
    if (params.runGeneAnnovar){ 
        ANNOVAR(
            input_ch, annodb, chr_prefix
        )
        logs     = logs.mix(ANNOVAR.out.log)
        versions = versions.mix(ANNOVAR.out.versions)
        input_ch = ANNOVAR.out.vcf
    }

    //
    // MODULE: SNV_RELIABILITY_PIPE
    //
    // RUN annotate_vcf.pl : BED files are used to annotate variants
    SNV_RELIABILITY_PIPE(
        input_ch, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
    )
    versions = versions.mix(SNV_RELIABILITY_PIPE.out.versions)

    // RUN: confidenceAnnotation_SNVs.py : Confidence annotation will be added to the variants
    vcf_ch = SNV_RELIABILITY_PIPE.out.vcf 
    input_ch = vcf_ch.join(ANNOTATE_VCF.out.median)
    
    CONFIDENCE_ANNOTATION(
        input_ch
    )
    ann_vcf_ch  = CONFIDENCE_ANNOTATION.out.vcf_ann
    versions    = versions.mix(CONFIDENCE_ANNOTATION.out.versions)

    //
    // MODULE: ANNOTATION_PIPES
    //
    // RUN annotate_vcf.pl : Uses optional databases to annotate variants, only given databases will be used. 
    //if (params.runSNVDeepAnnotation)
    //{
    //    ANNOTATION_PIPES (
    //    ann_vcf_ch, enchangers, cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets,
    //    cgi_mountains, phastconselem, encode_tfbs, mirnas_sncrnas
    //    )
    //    ann_vcf_ch  = ANNOTATION_PIPES.out.vcf 
    //    versions    = versions.mix(ANNOTATION_PIPES.out.versions)

    //}

    // ASK: If this is for the pancancer workflow, then also create a DKFZ specific file.// ask this
    // mkfifo is not implemented! If true runArtifactFilter creates a bias file will be used to plot errors
    // RUN REMANING AS ALL
    //if (params.runArtifactFilter){
    //        POST_PROCESS(
    //        ann_vcf_ch, ref 
    //    )
    //    versions = versions.mix(POST_PROCESS.out.versions)
    //}
    //else{
    //    FILTER_PEOVERLAP(
    //        ann_vcf_ch, ref  
    //    )
    //}

    //FILTER_PEOVERLAP(
    //    vcf_ch, ref
    //)
    //versions = versions.mix(FILTER_PEOVERLAP.out.versions)

    // createErrorPlots.py! works only if there is a bias file produced 
    //ERROR_PLOTS(
    //    FILTER_PEOVERLAP.out.somatic_snvs_for_bias 
    //)
    //versions = versions.mix(ERROR_PLOTS.out.versions)

    // plotBaseScoreDistribution.R That will work only if it is true
    //PLOT_BASESCORE_BIAS(
    //    ERROR_PLOTS.out.error_matrix
    //)


emit:
logs
versions
}
