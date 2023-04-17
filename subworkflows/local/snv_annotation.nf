//
// SNV ANNOTATION: USES A CUSTOM PERL SCRIPT TO ANNOTATE VCFS AND ANNOVAR
//

params.options = [:]

include { ANNOTATE_VCF             } from '../../modules/local/annotate_vcf.nf'          addParams( options: params.options )
include { ANNOVAR                  } from '../../modules/local/annovar.nf'               addParams( options: params.options )
include { SNV_RELIABILITY_PIPE     } from '../../modules/local/snv_reliability_pipe.nf'  addParams( options: params.options )
include { ANNOTATION_PIPES         } from '../../modules/local/annotation_pipes.nf'      addParams( options: params.options )
include { CONFIDENCE_ANNOTATION    } from '../../modules/local/confidence_annotation.nf' addParams( options: params.options )
include { FLAG_BIAS as FLAG_BIAS_1 } from '../../modules/local/flag_bias.nf'             addParams( options: params.options )
include { FLAG_BIAS as FLAG_BIAS_2 } from '../../modules/local/flag_bias.nf'             addParams( options: params.options )
include { TABIX_BGZIPTABIX         } from '../../modules/nf-core/modules/tabix/bgziptabix/main'            addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_1  } from '../../modules/local/filter_peoverlap.nf'         addParams( options: params.options )
include { FILTER_PEOVERLAP as FILTER_PEOVERLAP_2  } from '../../modules/local/filter_peoverlap.nf'         addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_1            } from '../../modules/local/error_plots.nf'              addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_2            } from '../../modules/local/error_plots.nf'              addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_3            } from '../../modules/local/error_plots.nf'              addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_4            } from '../../modules/local/error_plots.nf'              addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_1 } from '../../modules/local/plot_basescore_bias.nf' addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_2 } from '../../modules/local/plot_basescore_bias.nf' addParams( options: params.options )


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

    versions = Channel.empty()
    logs     = Channel.empty()
    plots_ch = Channel.empty() 

    //
    // MODULE: ANNOTATE_VCF
    //
    // RUN annotate_vcf.pl: Uses various databases (all mandatory exept recurrance) to annotate variants
    ANNOTATE_VCF (
        vcf_ch, 
        kgenome, 
        dbsnpsnv, 
        localcontrolwgs,
        localcontrolwes, 
        gnomadgenomes, 
        gnomadexomes, 
        chr_prefix
    )
    versions  = versions.mix(ANNOTATE_VCF.out.versions)
    // prepare input channel for annovar 
    // input_ch = meta, vcf, for_annovar.bed
    ch_vcf = ANNOTATE_VCF.out.unziped_vcf
    input_ch = ch_vcf.join(ANNOTATE_VCF.out.forannovar)

    //
    // MODULE: ANNOVAR
    // 
    // RUN annovar, processAnnovarOutput.pl and newCols2vcf.pl: annovar annotates and classifies the variants, 
    // perl scripts re-creates vcfs. 
    ANNOVAR(
        input_ch, 
        annodb, 
        chr_prefix
    )
    logs     = logs.mix(ANNOVAR.out.log)
    versions = versions.mix(ANNOVAR.out.versions)

    //
    // MODULE: SNV_RELIABILITY_PIPE
    //
    // RUN annotate_vcf.pl : BED files are used to annotate variants
    SNV_RELIABILITY_PIPE(
        ANNOVAR.out.vcf, 
        repeatmasker, 
        dacblacklist, 
        dukeexcluded, 
        hiseqdepth, 
        selfchain, 
        mapability, 
        simpletandemrepeats
    )
    versions = versions.mix(SNV_RELIABILITY_PIPE.out.versions)

    //
    // MODULE: CONFIDENCE_ANNOTATION
    //
    // Apply confidenceAnnotation_SNVs.py only to with germline variants
    CONFIDENCE_ANNOTATION(
        SNV_RELIABILITY_PIPE.out.vcf
    )
    versions = versions.mix(CONFIDENCE_ANNOTATION.out.versions)

    // If true runArtifactFilter creates a bias file will be used to plot errors
    if (params.runArtifactFilter){
        //
        // MODULE: FILTER_PEOVERLAP
        //
        FILTER_PEOVERLAP_1(
            CONFIDENCE_ANNOTATION.out.vcf, 
            ref
        )
        versions    = versions.mix(FILTER_PEOVERLAP_1.out.versions)
        altbasequal = FILTER_PEOVERLAP_1.out.alternative_allele_base_qualities
        refbasequal = FILTER_PEOVERLAP_1.out.reference_allele_base_qualities 
        altreadpos  = FILTER_PEOVERLAP_1.out.alternative_allele_read_positions 
        refreadpos  = FILTER_PEOVERLAP_1.out.reference_allele_read_positions

        // filter out the lists if no variant exists for visualization
        FILTER_PEOVERLAP_1.out.somatic_snvs
            .filter{meta, somatic_snvs -> WorkflowCommons.getNumLinesInFile(somatic_snvs) > 1}
            .set{somatic_vcf}

        ////////////////////////////////////////
        //// First round of plot generation ////
        ////////////////////////////////////////
        //
        // MODULE: ERROR_PLOTS
        //
        // Sequencing Error plot
        ERROR_PLOTS_1(
            somatic_vcf,
            'sequencing_specific', 
            'sequencing_specific_error_plot_before_filter', 
            'sequencing_error_matrix_first', 
            'Sequencing strand bias before guanine oxidation filter'
        )
        versions = versions.mix(ERROR_PLOTS_1.out.versions)
        plots_ch = plots_ch.mix(ERROR_PLOTS_1.out.plot)
        // Sequence Error plot
        ERROR_PLOTS_2(
            somatic_vcf, 
            'sequence_specific', 
            'sequence_specific_error_plot_before_filter',
            'sequence_error_matrix_first', 
            'PCR strand bias before guanine oxidation filter'
        )
        plots_ch = plots_ch.mix(ERROR_PLOTS_2.out.plot)

        //
        // MODULE: PLOT_BASESCORE_BIAS
        //
        // Run tripletBased_BQRatio_plotter.R only if generateExtendedQcPlots is true, this step only generates a pdf!
        if (params.generateExtendedQcPlots){
            // create input channel with error matrixes: 
            //Somatic SVC Temp,reference_allele_base_qualities, alternative_allele_base_qualities 
            somatic_vcf.join(refbasequal)
                .join(altbasequal)
                .set{ somatic_ch }
            
            PLOT_BASESCORE_BIAS_1(
                somatic_ch, 
                'base_score_bias_before_filter',
                'Base Quality Bias Plot for PID before guanine oxidation filter'
                )
            versions = versions.mix(PLOT_BASESCORE_BIAS_1.out.versions)
            plots_ch = plots_ch.mix(PLOT_BASESCORE_BIAS_1.out.plot) 
        }else{
            println "Extended QC plots not generated because generateExtendedQcPlots is set to ${params.generateExtendedQcPlots}"
        }

        //
        // MODULE: FLAG_BIAS
        //
        // create input channel for flag bias with error matrixes from error plots
        // error_ch: meta, _peoverlap.vcf, _sequence_error_matrix.txt, _sequencing_error_matrix.txt
        FILTER_PEOVERLAP_1.out.vcf.join(ERROR_PLOTS_2.out.error_matrix)
            .join(ERROR_PLOTS_1.out.error_matrix)
            .set{ error_ch }
        
        FLAG_BIAS_1(
            error_ch, 
            ref, 
            1
            )
        versions = versions.mix(FLAG_BIAS_1.out.versions)

        /////////////////////////////////////////    
        //// Second round of plot generation ////
        /////////////////////////////////////////    
        //
        // MODULE: ERROR_PLOTS
        //
        // Sequencing Error plot
        somatic_vcf = FLAG_BIAS_1.out.somatic_snvs
        ERROR_PLOTS_3(
            somatic_vcf,
            'sequencing_specific', 
            'sequencing_specific_error_plot_after_filter_once', 
            'sequencing_error_matrix_second', 
            'Sequencing strand bias after first round of guanine oxidation filter'
        )
        plots_ch = plots_ch.mix(ERROR_PLOTS_3.out.plot)

        // Sequence Error plot
        ERROR_PLOTS_4(
            somatic_vcf, 
            'sequence_specific', 
            'sequence_specific_error_plot_after_filter_once',
            'sequence_error_matrix_second', 
            'PCR strand bias after first round of guanine oxidation filter'
        )
        plots_ch = plots_ch.mix(ERROR_PLOTS_4.out.plot)

        //
        // MODULE: PLOT_BASESCORE_BIAS
        //
        // Run plot_basescore_bias.r only if generateExtendedQcPlots is true, this step only generates a pdf!
        if (params.generateExtendedQcPlots){
            // create input channel with error matrixes: somatic_vcf,reference_allele_base_qualities, alternative_allele_base_qualities 
            somatic_vcf.join(refbasequal)
                .join(altbasequal)
                .set{ som2_ch }

            PLOT_BASESCORE_BIAS_2(
                som2_ch, 
                'base_score_bias_after_filter_once',
                'Base Quality Bias Plot for PID after first round of guanine oxidation filter'
                )
            plots_ch = plots_ch.mix(PLOT_BASESCORE_BIAS_2.out.plot) 
        }else{
            println "Extended QC plots not generated because generateExtendedQcPlots is set to ${params.generateExtendedQcPlots}"
        }

        //
        // MODULE: FLAG_BIAS
        //
        // create input channel for flag bias with error matrixes from error plots
        FLAG_BIAS_1.out.vcf.join(ERROR_PLOTS_2.out.error_matrix)
            .join(ERROR_PLOTS_1.out.error_matrix)
            .set{ error2_ch }
        // input_ch: meta, _peoverlap.vcf, _sequence_error_matrix.txt, _sequencing_error_matrix.txt
        FLAG_BIAS_2(
            error2_ch, 
            ref, 
            2
            )        
        out_vcf = FLAG_BIAS_2.out.vcf

    }
    // IF runArticantfilter is false run only FILTER_PEOVERLAP
    else{
        println "QC filter not applied because runArticantfilter is set to ${params.runArticantfilter}"
        //
        // MODULE: FILTER_PEOVERLAP
        //
        FILTER_PEOVERLAP_2(
            CONFIDENCE_ANNOTATION.out.vcf, 
            ref 
        )
        out_vcf     = FILTER_PEOVERLAP_2.out.vcf
        altbasequal = FILTER_PEOVERLAP_2.out.alternative_allele_base_qualities
        refbasequal = FILTER_PEOVERLAP_2.out.reference_allele_base_qualities 
        altreadpos  = FILTER_PEOVERLAP_2.out.alternative_allele_read_positions 
        refreadpos  = FILTER_PEOVERLAP_2.out.reference_allele_read_positions     
    }
    
    //
    // MODULE: TABIX_BGZIPTABIX
    //
    TABIX_BGZIPTABIX(
        out_vcf
    )
    versions = versions.mix(TABIX_BGZIPTABIX.out.versions) 
    vcf_ch = TABIX_BGZIPTABIX.out.gz_tbi 

    //
    // MODULE: ANNOTATION_PIPES
    //
    // RUN annotate_vcf.pl : Uses optional databases to annotate variants, only given databases will be used. 
    if (params.runSNVDeepAnnotation)
    {
        ANNOTATION_PIPES (
            TABIX_BGZIPTABIX.out.gz_tbi, 
            enchangers, 
            cpgislands, 
            tfbscons, 
            encode_dnase, 
            mirnas_snornas, 
            cosmic, 
            mirbase, 
            mir_targets,
            cgi_mountains, 
            phastconselem, 
            encode_tfbs, 
            mirnas_sncrnas
        )
        vcf_ch   = ANNOTATION_PIPES.out.vcf 
        versions = versions.mix(ANNOTATION_PIPES.out.versions)
    }
    else{
        println "SNVDeep annotation not applied because runSNVDeepAnnotation is set to ${params.runSNVDeepAnnotation}"
    }

emit:
vcf_ch
altbasequal
refbasequal
altreadpos 
refreadpos
plots_ch
logs
versions
}
