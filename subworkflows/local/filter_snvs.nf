//
// FILTER SNVs: Filtering, extraction and plots
//

params.options = [:]

include { FILTER_BY_CRIT         } from '../../modules/local/filter_by_crit.nf'          addParams( options: params.options )
include { DBSNP_COUNTER          } from '../../modules/local/dbsnp_counter.nf'           addParams( options: params.options )
include { MUTATION_DISTANCE      } from '../../modules/local/mutation_distance.nf'       addParams( options: params.options )
include { INTERMUTATION_DISTANCE } from '../../modules/local/intermutation_distance.nf'  addParams( options: params.options )
include { PER_CHROM_PLOT         } from '../../modules/local/per_chrom_plot.nf'          addParams( options: params.options )
include { CONTEXT_FREQUENCIES    } from '../../modules/local/context_frequencies.nf'     addParams( options: params.options )
include { CONTEXT_PLOT           } from '../../modules/local/context_plot.nf'            addParams( options: params.options )
include { PURITY_RELOADED        } from '../../modules/local/purity_reloaded.nf'         addParams( options: params.options )
include { JSON_REPORT            } from '../../modules/local/json_report.nf'             addParams( options: params.options )
include { MERGE_PLOTS            } from '../../modules/local/merge_plots.nf'             addParams( options: params.options )
include { SNV_EXTRACTOR          } from '../../modules/local/snv_extractor.nf'           addParams( options: params.options )
include { TRIPLET_PLOTTER        } from '../../modules/local/triplet_plotter.nf'         addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_5     } from '../../modules/local/error_plots.nf'   addParams( options: params.options )
include { ERROR_PLOTS as ERROR_PLOTS_6     } from '../../modules/local/error_plots.nf'   addParams( options: params.options )
include { PLOT_BASESCORE_DISTRIBUTION      } from '../../modules/local/plot_basescore_distribution.nf'     addParams( options: params.options )
include { PLOT_BASESCORE_BIAS as PLOT_BASESCORE_BIAS_3 } from '../../modules/local/plot_basescore_bias.nf' addParams( options: params.options )
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_2       } from '../../modules/nf-core/modules/tabix/bgziptabix/main' addParams( options: params.options )

workflow FILTER_SNVS {
    take:
    input_ch        // channel: [val(meta), vcf, index,  altbasequal, refbasequal, altreadpos, refreadpos, plots ]
    ref             // reference channel [ref.fa, ref.fa.fai]
    chr_prefix      // val channel
    chrlength       // chr file    

    main:
    versions = Channel.empty()
    plots_ch = Channel.empty()

    // rawvcf_ch=meta, vcf, index
    rawvcf_ch  = input_ch.map{ it -> tuple( it[0], it[1], it[2] )}
    //
    // MODULE: FILTER_BY_CRIT
    //
    // RUN vcf_filter_bycrit.pl : filter only be apply on for no-control cases
    FILTER_BY_CRIT(
        rawvcf_ch
    )
    versions = versions.mix(FILTER_BY_CRIT.out.versions)

    //
    // MODULE: SNV_EXTRACTOR
    //
    // run snv_extractor_v1.pl
    SNV_EXTRACTOR(
        FILTER_BY_CRIT.out.vcf 
    )
    versions = versions.mix(SNV_EXTRACTOR.out.versions)

    // filter out the lists if no variant exists for visualization
    SNV_EXTRACTOR.out.somatic_snv
        .filter{meta, somatic_snv -> WorkflowCommons.getNumLinesInFile(somatic_snv) > 1}
        .set{orjinal_somatic_ch}

    somatic_ch = orjinal_somatic_ch
    // if rerun is false
    // Rest is the usual pipeline. if rerun is false
    if (params.runplots){

        // 1. Rainfall plots
        //
        // MODULE: DBSNP_COUNTER
        //
        // run  in_dbSNPcounter.pl     
        DBSNP_COUNTER(
            somatic_ch    
        )
        versions = versions.mix(DBSNP_COUNTER.out.versions)
        
        //
        // MODULE:JSON_REPORT
        //
        // prepare report json file
        // determine fraction of SNVs called as "synonymous SNV" among all exonic SNVs (QC value)
        // MAF plots are generated in the JSON_REPORT module
        // and THA detection
        JSON_REPORT(
            somatic_ch.join(DBSNP_COUNTER.out.indbsnp)
        )
        versions = versions.mix(JSON_REPORT.out.versions)
        plots_ch = plots_ch.mix(JSON_REPORT.out.plot)

        // Exclude Y is false. No gender info is collected. 
        //exy="--excludedChromosomes=chrY", ignoreY=1
        //mutation distance, rainfall plot,  mutation classes per chromosome
        //
        // MODULE: MUTATION_DISTANCE
        //
        // run mutationDistance.py
        MUTATION_DISTANCE(
            somatic_ch    
        )
        versions = versions.mix(MUTATION_DISTANCE.out.versions) 

        //
        // MODULE: INTERMUTATION_DISTANCE
        //
        // run intermutationDistance_Coord_color.r
        INTERMUTATION_DISTANCE(
            somatic_ch, 
            chrlength, 
            chr_prefix
        )
        versions = versions.mix(INTERMUTATION_DISTANCE.out.versions)
        plots_ch = plots_ch.mix(INTERMUTATION_DISTANCE.out.plot)

        //
        // MODULE: PER_CHROM_PLOT
        //   
        // run snvsPerChromPlot.r
        PER_CHROM_PLOT(
            MUTATION_DISTANCE.out.distance, 
            chrlength
        )  
        versions = versions.mix(PER_CHROM_PLOT.out.versions)
        plots_ch = plots_ch.mix(PER_CHROM_PLOT.out.plot) 

        //
        // MODULE: CONTEXT_FREQUENCIES
        //
        // Run SNV_context_frequencies.pl
        CONTEXT_FREQUENCIES(
            somatic_ch  
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
        plots_ch = plots_ch.mix(CONTEXT_PLOT.out.plot)

        //
        // MODULE: ERROR_PLOTS
        //
        // Sequencing Error plot
        ERROR_PLOTS_5(
            somatic_ch,
            'sequencing_specific', 
            "sequencing_specific_error_plot_conf_${params.min_confidence_score}_to_10", 
            "sequencing_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 
            'Final sequencing strand bias from vcf filter script'
        )
        plots_ch = plots_ch.mix(ERROR_PLOTS_5.out.plot)
        // Sequence Error plot
        ERROR_PLOTS_6(
            somatic_ch, 
            'sequence_specific',
            "sequence_specific_error_plot_conf_${params.min_confidence_score}_to_10",
            "sequence_specific_error_Matrix_conf_${params.min_confidence_score}_to_10", 
            'Final PCR strand bias from vcf filter script'
        )
        plots_ch = plots_ch.mix(ERROR_PLOTS_6.out.plot)
            
        if (params.generateExtendedQcPlots){
            //make base score bias and base score distribution plots
            //!!!!!! this step should only run if  refbasequal, altbasequal generated !!!!!!
            //
            // MODULE: PLOT_BASESCORE_BIAS
            //
            // run tripletBased_BQRatio_plotter.R
            // create input channel: meta, sometic_vcf,  refbasequal, altbasequal
            filt_ch  = input_ch.map{it -> tuple( it[0], it[4], it[3] )}

            PLOT_BASESCORE_BIAS_3(
                somatic_ch.join(filt_ch), 
                "base_score_bias_plot_conf_${params.min_confidence_score}_to_10", 
                "Final Base Quality Bias Plot for PID"
            )
            plots_ch = plots_ch.mix(PLOT_BASESCORE_BIAS_3.out.plot) 

            //
            // MODULE: PLOT_BASESCORE_DISTRIBUTION
            //
            //run plotBaseScoreDistribution.R
            PLOT_BASESCORE_DISTRIBUTION(
                somatic_ch.join(filt_ch),
                "base_score_distribution", 
                "for somatic SNVs for PID"
            )
            versions = versions.mix(PLOT_BASESCORE_DISTRIBUTION.out.versions)
            plots_ch = plots_ch.mix(PLOT_BASESCORE_DISTRIBUTION.out.plot)

            //
            // MODULE: TRIPLET_PLOTTER
            //
            // Run tripletBased_BQDistribution_plotter.R
            triplet_ch     = input_ch.map{ it -> tuple( it[0], it[3], it[4], it[5], it[6] )}
            
            TRIPLET_PLOTTER(
                somatic_ch.join(triplet_ch), 
                "Base score distribution of PID",
                0
            )
            versions = versions.mix(TRIPLET_PLOTTER.out.versions)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_1)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_2)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_3)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_4)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_5)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_6)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_7)
            plots_ch = plots_ch.mix(TRIPLET_PLOTTER.out.plot_8)
        }
        else{
            println "Extended QC plots not generated because generateExtendedQcPlots is set to ${params.generateExtendedQcPlots}"
        } 
    }
    else{
        println "Plots not generated because runplots is set to ${params.runplots}"
    }

    //
    // MODULE: MERGE_PLOTS
    // 
    // run ghostscript
    temp_ch = input_ch.map{ it -> tuple( it[0], it[7])} 
    temp2_ch = plots_ch.groupTuple().map {it -> tuple( it[0], it[1] )}            

    MERGE_PLOTS(
        temp_ch.join(temp2_ch)
    )
    versions = versions.mix(MERGE_PLOTS.out.versions) 

    // 3. Run purityEST
    if (params.runpurest){
        if (params.min_confidence_score > 7){
            //
            // MODULE: PURITY_RELOADED
            //
            // Run PurityReloaded.py
            PURITY_RELOADED(
                rawvcf_ch
            )
            versions = versions.mix(PURITY_RELOADED.out.versions) 
        }else{
            println "PurityEST not run because min_confidence_score is set to ${params.min_confidence_score} and purityEST requires a min_confidence_score of 8 or higher"
        }
    }
    
    emit:
    versions
}