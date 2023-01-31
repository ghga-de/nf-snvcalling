//
// FILTER SNVs: Filtering, extraction and plots
//

params.options = [:]

include { FILTER_BY_CRIT         } from '../../modules/local/filter_by_crit.nf'         addParams( options: params.options )
include { SNV_EXTRACTOR          } from '../../modules/local/snv_extractor.nf'          addParams( options: params.options )
include { DBSNP_COUNTER          } from '../../modules/local/dbsnp_counter.nf'          addParams( options: params.options )
include { MUTATION_DISTANCE      } from '../../modules/local/mutation_distance.nf'      addParams( options: params.options )
include { INTERMUTATION_DISTANCE } from '../../modules/local/intermutation_distance.nf' addParams( options: params.options )
include { PER_CHROM_PLOT         } from '../../modules/local/per_chrom_plot.nf'         addParams( options: params.options )
include { MAKE_MAF_INPUT         } from '../../modules/local/make_maf_input.nf'         addParams( options: params.options )
include { MAF_PLOTS              } from '../../modules/local/maf_plots.nf'              addParams( options: params.options )
include { CONTEXT_FREQUENCIES    } from '../../modules/local/context_frequencies.nf'    addParams( options: params.options )
include { CONTEXT_PLOT           } from '../../modules/local/context_plot.nf'           addParams( options: params.options )
include { PURITY_RELOADED        } from '../../modules/local/purity_reloaded.nf'        addParams( options: params.options )
include { PLOT_BASESCORE_DISTRIBUTION    } from '../../modules/local/plot_basescore_distribution.nf'           addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_5  } from '../../modules/local/error_plots.nf'     addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_6  } from '../../modules/local/error_plots.nf'     addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_3 } from '../../modules/local/plot_basescore_bias.nf'    addParams( options: params.options )

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

    //
    // MODULE: SNV_EXTRACTOR
    //

    SNV_EXTRACTOR(
    FILTER_BY_CRIT.out.vcf  
    )
    versions = versions.mix(SNV_EXTRACTOR.out.versions)

    // ASK about rerunning filter step

    if (params.rerunplots)
    {
        //!!!! track analysed chromosomes, if there is no Y, set ignoreY=1 and exY="--excludedChromosomes=chrY" 
        // 1. Rainfall plots
        //
        // MODULE: DBSNP_COUNTER
        //
        
        DBSNP_COUNTER(
        SNV_EXTRACTOR.out.somatic_snv    
        )
        versions = versions.mix(DBSNP_COUNTER.out.versions)  

        exy="--excludedChromosomes=chrY"
        //
        // MODULE: MUTATION_DISTANCE
        //

        MUTATION_DISTANCE(
        SNV_EXTRACTOR.out.somatic_snv, exy    
        )
        versions = versions.mix(MUTATION_DISTANCE.out.versions) 

        //
        // MODULE: INTERMUTATION_DISTANCE
        //
        // run intermutationDistance_Coord_color.r

        //
        // MODULE: PER_CHROM_PLOT
        //   
        // run snvsPerChromPlot.r

        // 2. MAF plots
        //
        // MODULE: MAKE_MAF_INPUT
        //
        // RUN makeMAFinput.pl
        MAKE_MAF_INPUT(
        SNV_EXTRACTOR.out.somatic_snv    
        )
        versions = versions.mix(MAKE_MAF_INPUT.out.versions)

        //
        // MODULE: MAF_PLOTS
        // 
        //Run MAF_plots.r
        //input_ch = meta, SNV_EXTRACTOR.out.somatic_snv, MAKE_MAF_INPUT.out.maf_values, DBSNP_COUNTER.out.indbsnp 
        input_ch = SNV_EXTRACTOR.out.somatic_snv.join(MAKE_MAF_INPUT.out.maf_values)
        input_ch = input_ch.join(DBSNP_COUNTER.out.indbsnp )
        MAF_PLOTS(
        input_ch
        )
        versions = versions.mix(MAF_PLOTS.out.versions)

        //
        // MODULE: CONTEXT_FREQUENCIES
        //
        // Run SNV_context_frequencies.pl
        CONTEXT_FREQUENCIES(
        SNV_EXTRACTOR.out.somatic_snv  
        )
        versions = versions.mix(CONTEXT_FREQUENCIES.out.versions)

        //
        // MODULE: CONTEXT_PLOT
        //
        // Run SNVSeqContext.R
        CONTEXT_PLOT(
        CONTEXT_FREQUENCIES.out.seq_contex  
        )
        versions = versions.mix(CONTEXT_PLOT.out.versions) 

        //
        // MODULE: ERROR_PLOTS
        //
        // Sequencing Error plot
        ERROR_PLOTS_5(
            SNV_EXTRACTOR.out.somatic_snv,'sequencing_specific', "sequencing_specific_error_plot_conf_${params.min_confidence_score}_to_10", "_sequencing_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 'Final sequencing strand bias from vcf filter script'
        )
        versions = versions.mix(ERROR_PLOTS_5.out.versions)

        // Sequence Error plot
        ERROR_PLOTS_6(
            SNV_EXTRACTOR.out.somatic_snv, 'sequence_specific',"sequence_specific_error_plot_conf_${params.min_confidence_score}_to_10","_sequence_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 'Final PCR strand bias from vcf filter script'
        )
        versions = versions.mix(ERROR_PLOTS_6.out.versions)

        if (params.generateExtendedQcPlots){

            //
            // MODULE: PLOT_BASESCORE_BIAS
            //

            //
            // MODULE: PLOT_BASESCORE_DISTIBUTION
            //
        }

    }
    // 3. Run purityEST
    if (params.runpurest){
        //
        // MODULE: PURITY_RELOADED
        //
        // Run PurityReloaded.py

        PURITY_RELOADED(
            vcf_ch
        )

    }

    emit:
    versions
}