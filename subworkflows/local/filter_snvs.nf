//
// FILTER SNVs: Filtering, extraction and plots
//

params.options = [:]

include { FILTER_BY_CRIT         } from '../../modules/local/filter_by_crit.nf'         addParams( options: params.options )
include { DBSNP_COUNTER          } from '../../modules/local/dbsnp_counter.nf'          addParams( options: params.options )
include { MUTATION_DISTANCE      } from '../../modules/local/mutation_distance.nf'      addParams( options: params.options )
include { INTERMUTATION_DISTANCE } from '../../modules/local/intermutation_distance.nf' addParams( options: params.options )
include { PER_CHROM_PLOT         } from '../../modules/local/per_chrom_plot.nf'         addParams( options: params.options )
include { MAKE_MAF_INPUT         } from '../../modules/local/make_maf_input.nf'         addParams( options: params.options )
include { MAF_PLOTS              } from '../../modules/local/maf_plots.nf'              addParams( options: params.options )
include { CONTEXT_FREQUENCIES    } from '../../modules/local/context_frequencies.nf'    addParams( options: params.options )
include { CONTEXT_PLOT           } from '../../modules/local/context_plot.nf'           addParams( options: params.options )
include { PURITY_RELOADED        } from '../../modules/local/purity_reloaded.nf'        addParams( options: params.options )
include { TRIPLET_PLOTTER        } from '../../modules/local/triplet_plotter.nf'        addParams( options: params.options )
include { BEDTOOLS_SUBTRACT      } from '../../modules/local/bedtools_subtract.nf'      addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_5     } from '../../modules/local/error_plots.nf'   addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_6     } from '../../modules/local/error_plots.nf'   addParams( options: params.options )
include { SNV_EXTRACTOR as SNV_EXTRACTOR_1 } from '../../modules/local/snv_extractor.nf' addParams( options: params.options )
include { SNV_EXTRACTOR as SNV_EXTRACTOR_2 } from '../../modules/local/snv_extractor.nf' addParams( options: params.options )
include { PLOT_BASESCORE_DISTRIBUTION      } from '../../modules/local/plot_basescore_distribution.nf'     addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_3 } from '../../modules/local/plot_basescore_bias.nf' addParams( options: params.options )

workflow FILTER_SNVS {
    take:
    input_ch        // channel: [val(meta), vcf, index,  altbasequal, refbasequal, altreadpos, refreadpos ]
    ref             // reference channel [ref.fa, ref.fa.fai]
    chr_prefix      // val channel
    chrlength       // chr file    

    main:

    versions=Channel.empty()

    //
    // MODULE: FILTER_BY_CRIT
    //
    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
    input_ch
    )
    versions = versions.mix(FILTER_BY_CRIT.out.versions)

    //
    // MODULE: SNV_EXTRACTOR
    //

    SNV_EXTRACTOR_1(
    FILTER_BY_CRIT.out.vcf  
    )
    somatic_vcf_ch = SNV_EXTRACTOR_1.out.somatic_snv
    versions = versions.mix(SNV_EXTRACTOR_1.out.versions)

    // Rerun Filtering
    if (params.rerunfiltering)
    {
        ///!!! This part is not working for now !!!!!

        // temp_ch=meta, somatic_vcf, vcfgz, index,  altbasequal, refbasequal, altreadpos, refreadpos
        temp_ch =somatic_vcf_ch.join(input_ch)

        //
        // MODULE: TRIPLET_PLOTTER
        //
        // Run tripletBased_BQDistribution_plotter.R

        TRIPLET_PLOTTER(
            temp_ch, 1, "Base score distribution of PID \nafter Median'${params.median_filter_threshold}' filtering"
        )
        versions = version.mix(TRIPLET_PLOTTER.out.versions)

        //
        // MODULE: SNV_EXTRACTOR
        //
        SNV_EXTRACTOR_2(
        TRIPLET_PLOTTER.out.filtered_vcf 
        )
        versions = versions.mix(SNV_EXTRACTOR_2.out.versions)

        //
        // MODULE: BEDTOOLS_SUBTRACT
        //

        temp3_ch = SNV_EXTRACTOR_2.out.somatic_snv.join(TRIPLET_PLOTTER.out.filtered_vcf)
        BEDTOOLS_SUBTRACT(
            temp3_ch
        )
        versions = versions.mix(BEDTOOLS_SUBTRACT.out.versions)

    }
    else {
        // Rest is the usual pipeline. if rerun is false
            if (params.runplots)
            {
            //!!!! track analysed chromosomes, if there is no Y, set ignoreY=1 and exY="--excludedChromosomes=chrY" 
            // 1. Rainfall plots
            //
            // MODULE: DBSNP_COUNTER
            //
        
            DBSNP_COUNTER(
            somatic_vcf_ch    
            )
            versions = versions.mix(DBSNP_COUNTER.out.versions)  

            exy="--excludedChromosomes=chrY"
            //
            // MODULE: MUTATION_DISTANCE
            //

            MUTATION_DISTANCE(
            somatic_vcf_ch, exy    
            )
            versions = versions.mix(MUTATION_DISTANCE.out.versions) 

            //
            // MODULE: INTERMUTATION_DISTANCE
            //
            // run intermutationDistance_Coord_color.r
            INTERMUTATION_DISTANCE(
                somatic_vcf_ch, chrlength, chr_prefix
            )
            versions = versions.mix(INTERMUTATION_DISTANCE.out.versions)
            plot1 = INTERMUTATION_DISTANCE.out.plot

            //
            // MODULE: PER_CHROM_PLOT
            //   
            // run snvsPerChromPlot.r
            PER_CHROM_PLOT(
                MUTATION_DISTANCE.out.distance, chrlength
            )  
            versions = versions.mix(PER_CHROM_PLOT.out.versions)
            plot2 = PER_CHROM_PLOT.out.plot         

            // 2. MAF plots
            //
            // MODULE: MAKE_MAF_INPUT
            //
            // RUN makeMAFinput.pl
            MAKE_MAF_INPUT(
            somatic_vcf_ch    
            )
            versions = versions.mix(MAKE_MAF_INPUT.out.versions)

            //
            // MODULE: MAF_PLOTS
            // 
            //Run MAF_plots.r
            //temp2_ch = meta, SNV_EXTRACTOR.out.somatic_snv, MAKE_MAF_INPUT.out.maf_values, DBSNP_COUNTER.out.indbsnp 
            temp2_ch = somatic_vcf_ch.join(MAKE_MAF_INPUT.out.maf_values)
            temp2_ch = temp2_ch.join(DBSNP_COUNTER.out.indbsnp )
            MAF_PLOTS(
            temp2_ch
            )
            versions = versions.mix(MAF_PLOTS.out.versions)
            plot3 = MAF_PLOTS.out.plot

            //
            // MODULE: CONTEXT_FREQUENCIES
            //
            // Run SNV_context_frequencies.pl
            CONTEXT_FREQUENCIES(
            somatic_vcf_ch  
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
            plot4 = CONTEXT_PLOT.out.plot

            //
            // MODULE: ERROR_PLOTS
            //
            // Sequencing Error plot
            ERROR_PLOTS_5(
            somatic_vcf_ch,'sequencing_specific', "sequencing_specific_error_plot_conf_${params.min_confidence_score}_to_10", "_sequencing_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 'Final sequencing strand bias from vcf filter script'
            )
            versions = versions.mix(ERROR_PLOTS_5.out.versions)
            plot5 = ERROR_PLOTS_5.out.plot

            // Sequence Error plot
            ERROR_PLOTS_6(
            somatic_vcf_ch, 'sequence_specific',"sequence_specific_error_plot_conf_${params.min_confidence_score}_to_10","_sequence_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 'Final PCR strand bias from vcf filter script'
            )
            versions = versions.mix(ERROR_PLOTS_6.out.versions)
            plot6 = ERROR_PLOTS_6.out.plot

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
            input_ch
            )
        // final json report is missing

        }
    }

    emit:
    versions
}