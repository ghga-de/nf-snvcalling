/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnvcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, 
                         params.fasta,
                         params.multiqc_config]

def checkPathParamList_annotation = [params.k_genome,
                                    params.dbsnp_snv,
                                    params.mapability_file,
                                    params.repeat_masker,
                                    params.simple_tandemrepeats]

for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// for annotate_vcf.pl 
if (params.runSNVAnnotation){
    for (param in checkPathParamList_annotation) { if (param) { file(param, checkIfExists: true) } }
}

// If runIndelDeepAnnotation is true; at least one of the annotation files must be provided
if ((params.runSNVDeepAnnotation) && (!params.enchancer_file && !params.cpgislands_file && !params.tfbscons_file && !params.encode_dnase_file && !params.mirnas_snornas_file && !params.mirna_sncrnas_file && !params.mirbase_file && !params.cosmic_file && !params.mir_targets_file && !params.cgi_mountains_file && !params.phastconselem_file && !params.encode_tfbs_file)) { 
    log.error "Please specify at least one annotation file to perform SNV Deep Annotation"
    exit 1
}

//
// Check mandatory parameters
//

if (params.input)         { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// Annovar only be checked if annovar is true
if (params.annotation_tool.contains("annovar")){
    file(params.annovar_path, checkIfExists: true)
}

if (params.annotation_tool.contains("vep")){
    if(!params.vep_cache_version || !params.vep_genome ){
        log.error "Please specify params.vep_cache_version and params.vep_genome to run VEP tool"
        exit 1
    }
    if (params.download_cache){
        file(params.vep_cache, checkIfExists: true)   
    }
}


// Set up reference depending on the genome choice
ref            = Channel.fromPath([params.fasta,params.fasta_fai], checkIfExists: true).collect()
chr_prefix     = params.chr_prefix  ? Channel.value(params.chr_prefix) : Channel.value("")
chrlength      = params.chrom_sizes ? Channel.fromPath(params.chrom_sizes, checkIfExists: true) : Channel.empty()   
contigs        = params.contig_file ? Channel.fromPath(params.contig_file, checkIfExists: true).map{it -> ["contigs", it]} : Channel.value([[],[]])
config         = Channel.fromPath("${projectDir}/assets/config/convertToStdVCF.json", checkIfExists: true).collect()

// Annovar table folder

annodb              = params.annovar_path       ? Channel.fromPath(params.annovar_path + '/humandb/') 
                                                : Channel.empty()
// VEP cache
vep_cache_db        = params.vep_cache          ? Channel.fromPath(params.vep_cache).collect()         : []

// Annotation databases
kgenome             = params.k_genome           ? Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
dbsnpsnv            = params.dbsnp_snv          ? Channel.fromPath([params.dbsnp_snv, params.dbsnp_snv + '.tbi'], checkIfExists: true).collect()  
                                                : Channel.value([[],[]])                                       
localcontrolwgs     = params.local_control_wgs  ? Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect()  
                                                : Channel.value([[],[]])
localcontrolwes     = params.local_control_wes  ? Channel.fromPath([params.local_control_wes,params.local_control_wes + '.tbi' ], checkIfExists: true).collect()     
                                                : Channel.value([[],[]])
gnomadgenomes       = params.gnomad_genomes     ? Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect() 
                                                : Channel.value([[],[]])
gnomadexomes        = params.gnomad_exomes      ? Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect()     
                                                : Channel.value([[],[]])


// Realiability files
repeatmasker        = params.repeat_masker      ? Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
dacblacklist        = params.dac_blacklist      ? Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
dukeexcluded        = params.duke_excluded      ? Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
hiseqdepth          = params.hiseq_depth        ? Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
selfchain           = params.self_chain         ? Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
mapability          = params.mapability_file    ? Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
simpletandemrepeats = params.simple_tandemrepeats ? Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
        
// Indel Deep Annotation files
enchangers          = params.enchancer_file     ? Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
cpgislands          = params.cpgislands_file    ? Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
tfbscons            = params.tfbscons_file      ? Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
encode_dnase        = params.encode_dnase_file  ? Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect() 
                                                : Channel.value([[],[]])
mirnas_snornas      = params.mirnas_snornas_file  ? Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
mirnas_sncrnas      = params.mirna_sncrnas_file ? Channel.fromPath([params.mirna_sncrnas_file, params.mirna_sncrnas_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
cosmic              = params.cosmic_file        ? Channel.fromPath([params.cosmic_file, params.cosmic_file + '.tbi'], checkIfExists: true).collect() 
                                                : Channel.value([[],[]])
mirbase             = params.mirbase_file       ? Channel.fromPath([params.mirbase_file, params.mirbase_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
mir_targets         = params.mir_targets_file   ? Channel.fromPath([params.mir_targets_file, params.mir_targets_file + '.tbi'], checkIfExists: true).collect() 
                                                : Channel.value([[],[]])
cgi_mountains       = params.cgi_mountains_file ? Channel.fromPath([params.cgi_mountains_file, params.cgi_mountains_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
phastconselem       = params.phastconselem_file ? Channel.fromPath([params.phastconselem_file, params.phastconselem_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
encode_tfbs         = params.encode_tfbs_file   ? Channel.fromPath([params.encode_tfbs_file, params.encode_tfbs_file + '.tbi'], checkIfExists: true).collect()
                                                : Channel.value([[],[]])
                                    
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK         } from '../subworkflows/local/input_check'
include { MPILEUP_SNV_CALL    } from '../subworkflows/local/mpileup_snv_call'
include { SNV_ANNOTATION      } from '../subworkflows/local/snv_annotation'
include { FILTER_SNVS         } from '../subworkflows/local/filter_snvs'
include { OUTPUT_STANDARD_VCF } from '../subworkflows/local/output_standard_vcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

//
// MODULE: Local Modules
//

include { GREP_SAMPLENAME   } from '../modules/local/grep_samplename.nf'
include { GETCHROMSIZES     } from '../modules/local/getchromsizes.nf'
include { GET_CONTIGS       } from '../modules/local/get_contigs.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SNVCALLING {

    ch_versions = Channel.empty()
    ch_logs     = Channel.empty()
    ch_stdvcf   = Channel.empty()
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    sample_ch   = INPUT_CHECK.out.ch_sample

    if ( !params.chrom_sizes) {
        //
        // MODULE: Prepare chromosome size file if not provided
        //
        GETCHROMSIZES(
            ref
            )
        ch_versions = ch_versions.mix(GETCHROMSIZES.out.versions)
        chrlength   = GETCHROMSIZES.out.sizes
    }
    interval_ch  = chrlength.splitCsv(sep: '\t', by:1)

    if (params.runcontigs != "NONE") {
        //
        // MODULE: Prepare contigs file if not provided
        //
        GET_CONTIGS(
            sample_ch,
            contigs
            )
        ch_versions = ch_versions.mix(GET_CONTIGS.out.versions)
        GET_CONTIGS.out.contigs.filter{meta, contig -> WorkflowCommons.getNumLinesInFile(contig) > 0}
                .set{contigs}
    }
    //
    // MODULE: Extract sample name from BAM
    //
    GREP_SAMPLENAME(
        sample_ch
    )
    ch_versions = ch_versions.mix(GREP_SAMPLENAME.out.versions)

    // Prepare an input channel of sample with sample names
    name_ch   = GREP_SAMPLENAME.out.samplenames
    ch_sample = sample_ch.join(name_ch)


    //
    // SUBWORKFLOW: MPILEUP_SNV_CALL: Call SNVs
    //
    MPILEUP_SNV_CALL(
        ch_sample, 
        ref, 
        interval_ch,
        contigs 
    )
    ch_versions = ch_versions.mix(MPILEUP_SNV_CALL.out.versions)

    //
    // SUBWORKFLOW: SNV_ANNOTATION: Annotate SNVs
    //
    
    if (params.runSNVAnnotation){ 
        SNV_ANNOTATION(
            MPILEUP_SNV_CALL.out.vcf_ch, 
            ref, 
            kgenome,dbsnpsnv,localcontrolwgs,localcontrolwes,gnomadgenomes,gnomadexomes,
            repeatmasker, dacblacklist, dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats,
            enchangers, cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets, cgi_mountains, phastconselem, encode_tfbs, mirnas_sncrnas, 
            chr_prefix,
            annodb,
            vep_cache_db
        )
        ch_versions = ch_versions.mix(SNV_ANNOTATION.out.versions)

        vcf_ch     = SNV_ANNOTATION.out.vcf_ch
        ch_stdvcf  = ch_stdvcf.mix(vcf_ch.map{ it -> tuple( it[0], it[1] )})
        

        //
        // SUBWORKFLOW: FILTER_SNVS: Filters SNVs
        //
        // input_ch= meta, annotated vcf, index, altbasequal, refbasequal, altreadpos, refreadpos, 
                    //sequence_spesific_error_plot_1, sequencing_spesific_error_plot_1, sequence_spesific_error_plot_2
                    //sequencing_spesific_error_plot_2, base_score_distribution_plot_1, base_score_distribution_plot_2
        vcf_ch.join(SNV_ANNOTATION.out.altbasequal)
                .join(SNV_ANNOTATION.out.refbasequal)
                .join(SNV_ANNOTATION.out.altreadpos)
                .join(SNV_ANNOTATION.out.refreadpos)
                .join(SNV_ANNOTATION.out.plots_ch.groupTuple())
                .set{input_ch}
        
        if (params.runSNVVCFFilter){
            FILTER_SNVS(
                input_ch, 
                ref, 
                chr_prefix, 
                chrlength
            )
            ch_versions = ch_versions.mix(FILTER_SNVS.out.versions)
            ch_stdvcf  = ch_stdvcf.mix(FILTER_SNVS.out.convert_snvs)
        }
        else{
            println "Skipping SNV filtering"
        }
    }
    else{
        println "Skipping SNV annotation and filtering"
    }


    if (params.standard_vcf){
        println "VCF output is standardizing.."

        OUTPUT_STANDARD_VCF(
            ch_stdvcf.combine(MPILEUP_SNV_CALL.out.vcf_ch, by:0),
            config
        )
        ch_versions = ch_versions.mix(OUTPUT_STANDARD_VCF.out.versions)
    }   

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    if (!params.skip_multiqc){
        //
        // MODULE: MultiQC
        //
        workflow_summary    = WorkflowSnvcalling.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        ch_multiqc_files = Channel.empty()
        ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
        ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        MULTIQC (
            ch_multiqc_files.collect()
        )
        multiqc_report = MULTIQC.out.report.toList()
        ch_versions    = ch_versions.mix(MULTIQC.out.versions)
    }
    else{
        println "Skipping MultiQC"
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
