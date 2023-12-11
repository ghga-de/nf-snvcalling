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
include { ENSEMBLVEP_VEP         } from '../../modules/nf-core/modules/ensemblvep/vep/main'       addParams( options: params.options )
include { ENSEMBLVEP_DOWNLOAD    } from '../../modules/nf-core/modules/ensemblvep/download/main'  addParams( options: params.options )


workflow SNV_ANNOTATION {
    take:
    vcf_ch           // channel: [val(meta), vcf.gz, vcf.gz.tbi  ]
    ref              // channel: [path(fasta), path(fai)]
    annotate_vcf_ref // channel: [val(meta2),file(kgenome),file(kgenome_i),file(dbsnpsnv),file(dbsnpsnv_i),file(localcontrolwgs),file(localcontrolwgs_i),file(localcontrolwes),file(localcontrolwes_i),file(gnomadgenomes),file(gnomadgenomes_i),file(gnomadexomes),file(gnomadexomes_i)]
    realibility_ref  // channel: [val(meta2),file(repeatmasker),file(repeatmasker_i),file(dacblacklist),file(dacblacklist_i),file(dukeexcluded),file(dukeexcluded_i),file(hiseqdepth),file(hiseqdepth_i),file(selfchain),file(selfchain_i),file(mapability),file(mapability_i),file(simpletandemrepeats),file(simpletandemrepeats_i)]
    deepanno_ref     // channel: [val(meta2),file(enchangers),file(enchangers_i),file(cpgislands),file(cpgislands_i),file(tfbscons),file(tfbscons_i),tuple file(encode_dnase),file(encode_dnase_i),file(mirnas_snornas),file(mirnas_snornas_i),file(cosmic),file(cosmic_i),file(mirbase),file(mirbase_i),file(mir_targets),file(mir_targets_i),file(cgi_mountains),file(cgi_mountains_i),file(phastconselem),file(phastconselem_i),file(encode_tfbs),file(encode_tfbs_i),file(mirnas_sncrnas),file(mirnas_sncrnas_i)]
    chr_prefix       // val channel: [prefix]
    annodb           // path: annovar db
    vep_cache        // path: vep cache

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
        annotate_vcf_ref, 
        chr_prefix
    )
    versions  = versions.mix(ANNOTATE_VCF.out.versions)

    ANNOTATE_VCF.out.unziped_vcf
        .join(ANNOTATE_VCF.out.forannovar)
        .set{anno_ch}

    if (params.annotation_tool.contains("annovar")){
        //
        // MODULE: ANNOVAR
        // 
        // RUN annovar, processAnnovarOutput.pl and newCols2vcf.pl: annovar annotates and classifies the variants, 
        // perl scripts re-creates vcfs. 
        ANNOVAR(
            anno_ch, 
            annodb, 
            chr_prefix
        )
        logs     = logs.mix(ANNOVAR.out.log)
        versions = versions.mix(ANNOVAR.out.versions)
        annotated_vcf = ANNOVAR.out.vcf
    }
    else{
        if(params.download_cache){
            ENSEMBLVEP_DOWNLOAD(
                anno_ch.map{ it -> tuple( it[0], it[1])}
                )
            versions  = versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)
            vep_cache = ENSEMBLVEP_DOWNLOAD.out.cache
        }
        
        ENSEMBLVEP_VEP(
            ANNOTATE_VCF.out.unziped_vcf,
            vep_cache,
            ref
        )
        versions      = versions.mix(ENSEMBLVEP_VEP.out.versions)
        annotated_vcf = ENSEMBLVEP_VEP.out.vcf
    }

    //
    // MODULE: SNV_RELIABILITY_PIPE
    //
    // RUN annotate_vcf.pl : BED files are used to annotate variants
    SNV_RELIABILITY_PIPE(
        annotated_vcf, 
        realibility_ref
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
            // Somatic SVC Temp,reference_allele_base_qualities, alternative_allele_base_qualities 
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
            deepanno_ref
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
