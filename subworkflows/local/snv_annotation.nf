//
// SNV ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF           } from '../../modules/local/annotate_vcf.nf'            addParams( options: params.options )
include { ANNOVAR                } from '../../modules/local/annovar.nf'                 addParams( options: params.options )
include { SNV_RELIABILITY_PIPE   } from '../../modules/local/snv_reliability_pipe.nf'    addParams( options: params.options )
include { ANNOTATION_PIPES       } from '../../modules/local/annotation_pipes.nf'        addParams( options: params.options )
include { CONFIDENCE_ANNOTATION} from '../../modules/local/confidence_annotation_1.nf' addParams( options: params.options )
include { TABIX_BGZIPTABIX       } from '../../modules/nf-core/modules/tabix/bgziptabix/main' addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_1  } from '../../modules/local/filter_peoverlap.nf'        addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_2  } from '../../modules/local/filter_peoverlap.nf'        addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_3  } from '../../modules/local/filter_peoverlap.nf'        addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_4  } from '../../modules/local/filter_peoverlap.nf'        addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_1            } from '../../modules/local/error_plots.nf'             addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_2            } from '../../modules/local/error_plots.nf'             addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_3            } from '../../modules/local/error_plots.nf'             addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_4            } from '../../modules/local/error_plots.nf'             addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_1 } from '../../modules/local/plot_basescore_bias.nf'    addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_2 } from '../../modules/local/plot_basescore_bias.nf'    addParams( options: params.options )
include { FLAG_BIAS as FLAG_BIAS_1 }        from '../../modules/local/flag_bias.nf'    addParams( options: params.options )
include { FLAG_BIAS as FLAG_BIAS_2 }        from '../../modules/local/flag_bias.nf'    addParams( options: params.options )

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
    // this is not correct!!!!!!!! should be fixed
    //
    // MODULE: ANNOVAR
    //
    if (params.runGeneAnnovar){ 
        ANNOVAR(
            input_ch, annodb, chr_prefix
        )
        logs     = logs.mix(ANNOVAR.out.log)
        versions = versions.mix(ANNOVAR.out.versions)
        ch_vcf   = ANNOVAR.out.vcf
        input_ch = ch_vcf.join(ANNOTATE_VCF.out.forannovar) 
    }

    //
    // MODULE: SNV_RELIABILITY_PIPE
    //
    // RUN annotate_vcf.pl : BED files are used to annotate variants
    SNV_RELIABILITY_PIPE(
        input_ch, repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats
    )
    versions = versions.mix(SNV_RELIABILITY_PIPE.out.versions)

    //
    // MODULE: CONFIDENCE_ANNOTATION_1
    //
    CONFIDENCE_ANNOTATION(
        SNV_RELIABILITY_PIPE.out.vcf
    )

    versions = versions.mix(CONFIDENCE_ANNOTATION.out.versions)

    // ASK: If this is for the pancancer workflow, then also create a DKFZ specific file.// ask this
    // mkfifo is not implemented! If true runArtifactFilter creates a bias file will be used to plot errors
    if (params.runArtifactFilter){
        FILTER_PEOVERLAP_1(
            CONFIDENCE_ANNOTATION.out.vcf, ref 
        )
        versions = versions.mix(FILTER_PEOVERLAP_1.out.versions)

        ERROR_PLOTS_1(
            FILTER_PEOVERLAP_1.out.somatic_snvs_tmp,'sequencing_specific', 'sequencing_specific_error_plot_before_filter', 'sequencing_error_matrix', 'Sequencing strand bias before guanine oxidation filter'
        )
        versions = versions.mix(ERROR_PLOTS_1.out.versions)

        ERROR_PLOTS_2(
            FILTER_PEOVERLAP_1.out.somatic_snvs_tmp, 'sequence_specific', 'sequence_specific_error_plot_before_filter','sequence_error_matrix', 'PCR strand bias before guanine oxidation filter'
        )
        versions = versions.mix(ERROR_PLOTS_2.out.versions)

        //
        // MODULE: PLOT_BASESCORE_BIAS
        //
        // Run plot_basescore_bias.r only if generateExtendedQcPlots is true, this step only generates a pdf!
        if (params.generateExtendedQcPlots){
            // create input channel with error matrixes
            input_ch = FILTER_PEOVERLAP_1.out.somatic_snvs_tmp.join(FILTER_PEOVERLAP_1.out.reference_allele_base_qualities)
            input_ch = input_ch.join(FILTER_PEOVERLAP_1.out.alternative_allele_base_qualities)

            PLOT_BASESCORE_BIAS_1(
                input_ch, 'base_score_bias_before_filter','Base Quality Bias Plot for PID before guanine oxidation filter'
                )
            versions = versions.mix(PLOT_BASESCORE_BIAS_1.out.versions)
        }

        // create input channel for flag bias with error matrixes from error plots
        input_ch = FILTER_PEOVERLAP_1.out.vcf.join(ERROR_PLOTS_2.out.error_matrix)
        input_ch = input_ch.join(ERROR_PLOTS_1.out.error_matrix)

        FLAG_BIAS_1(
            input_ch, ref
            )
        
    }
    // IF runArticantfilter is false run only FILTER_PEOVERLAP
    else{
        FILTER_PEOVERLAP_2(
            CONFIDENCE_ANNOTATION.out.vcf, ref  
        )
        versions = versions.mix(FILTER_PEOVERLAP_2.out.versions)

        TABIX_BGZIPTABIX(
            FILTER_PEOVERLAP_2.out.vcf
            )
        versions = versions.mix(TABIX_BGZIPTABIX.out.versions)       
    }

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


emit:
logs
versions
}
