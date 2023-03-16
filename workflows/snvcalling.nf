/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnvcalling.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config]

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
    log.error "Please specify at least one annotation file to perform INDEL Deep Annotation"
    exit 1
}

//
// Check mandatory parameters
//

if (params.input)         { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
// Annovar only be checked if runGeneAnnovar is true
if ((params.runCytoband ) &&(params.runGeneAnnovar) &&(params.annovar_path))  
    { annodb = Channel.fromPath(params.annovar_path + '/humandb/', checkIfExists: true ) } 
else 
    { annodb = Channel.empty() }

// Set up reference depending on the genome choice
// NOTE: link will be defined by aoutomatic reference generation when the pipeline ready!
if ((params.reference) && (params.chrlength) && (params.chr_prefix))
    {
        fa_file    = params.reference
        chr_file   = params.chrlength
        chr_prefix = params.chr_prefix

        if (params.contig_file) {
            contig_file = params.contig_file
        }
    }
    else{
    if (params.ref_type == 'hg37')
        { 
        fa_file  = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef_Phix/hs37d5_PhiX.fa"
        chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hs37d5.fa.chrLenOnlyACGT_realChromosomes.tab'      
        chr_prefix   = Channel.value("")                    
        }
    if (params.ref_type == 'hg19') 
        { 
        fa_file  = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_1-22_X_Y_M.fa"
        chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_1-22_X_Y_M.fa.chrLenOnlyACGT.tab'
        chr_prefix   = Channel.value("chr")
        }
    if (params.ref_type == 'hg38') 
        { 
        fa_file = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg_GRCh38/sequence/GRCh38_decoy_ebv_alt_hla_phiX.fa"
        chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg_GRCh38/stats/GRCh38_decoy_ebv_alt_hla_phiX.fa.chrLenOnlyACGT_realChromosomes.tsv'
        chr_prefix   = Channel.value("chr")
        // HLA and ALT contigs
        contig_file = 'assets/GRCh38_decoy_ebv_phiX_alt_hla_chr.fa.contig.bed'      
        }
    }

// prepare  channels
ref          = Channel.fromPath([fa_file,fa_file +'.fai'], checkIfExists: true).collect()
chrlength    = Channel.fromPath(chr_file, checkIfExists: true)
interval_ch  = chrlength.splitCsv(sep: '\t', by:1)
if (params.ref_type == 'hg38') { contigs = Channel.fromPath(contig_file, checkIfExists: true) } else { contigs = Channel.empty() }

// Annotation databases
if (params.k_genome)             { kgenome = Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect() } else { kgenome = Channel.of([],[]) }
if (params.dbsnp_snv)            { dbsnpsnv = Channel.fromPath([params.dbsnp_snv,params.dbsnp_snv +'.tbi' ], checkIfExists: true).collect() } else { dbsnpsnv = Channel.of([],[]) }
if (params.local_control_wgs)    { localcontrolwgs = Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwgs = Channel.of([],[]) }
if (params.local_control_wes)    { localcontrolwes = Channel.fromPath([params.local_control_wes,params.local_control_wes + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwes = Channel.of([],[]) }
if (params.gnomad_genomes)       { gnomadgenomes = Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect() } else { gnomadgenomes = Channel.of([],[]) }
if (params.gnomad_exomes)        { gnomadexomes = Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect() } else { gnomadexomes = Channel.of([],[]) }
// Annovar table folder
if (params.annovar_path)         { annodb = Channel.fromPath(params.annovar_path + '/humandb/', checkIfExists: true ) } else { annodb = Channel.empty() }
// Realiability files
if (params.repeat_masker)        { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.of([],[]) }
if (params.dac_blacklist)        { dacblacklist = Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect() } else { dacblacklist = Channel.of([],[]) }
if (params.duke_excluded)        { dukeexcluded = Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect() } else { dukeexcluded = Channel.of([],[]) }
if (params.hiseq_depth)          { hiseqdepth = Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect() } else { hiseqdepth = Channel.of([],[]) }
if (params.self_chain)           { selfchain = Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect() } else { selfchain = Channel.of([],[]) }
if (params.mapability_file)      { mapability = Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect() } else { mapability = Channel.of([],[]) }
if (params.simple_tandemrepeats) { simpletandemrepeats = Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect() } else { simpletandemrepeats = Channel.of([],[]) }
// Indel Deep Annotation files
if (params.enchancer_file)       { enchangers = Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect() } else { enchangers = Channel.of([],[]) }
if (params.cpgislands_file)      { cpgislands = Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect() } else { cpgislands = Channel.of([],[]) }
if (params.tfbscons_file)        { tfbscons = Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect() } else { tfbscons = Channel.of([],[]) }
if (params.encode_dnase_file)    { encode_dnase = Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect() } else { encode_dnase = Channel.of([],[]) }
if (params.mirnas_snornas_file)  { mirnas_snornas = Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect() } else { mirnas_snornas = Channel.of([],[]) }
if (params.mirna_sncrnas_file)   { mirnas_sncrnas = Channel.fromPath([params.mirna_sncrnas_file, params.mirna_sncrnas_file + '.tbi'], checkIfExists: true).collect() } else { mirnas_sncrnas = Channel.of([],[]) }
if (params.cosmic_file)          { cosmic = Channel.fromPath([params.cosmic_file, params.cosmic_file + '.tbi'], checkIfExists: true).collect() } else { cosmic = Channel.of([],[]) }
if (params.mirbase_file)         { mirbase = Channel.fromPath([params.mirbase_file, params.mirbase_file + '.tbi'], checkIfExists: true).collect() } else { mirbase = Channel.of([],[]) }
if (params.mir_targets_file)     { mir_targets = Channel.fromPath([params.mir_targets_file, params.mir_targets_file + '.tbi'], checkIfExists: true).collect() } else { mir_targets = Channel.of([],[]) }
if (params.cgi_mountains_file)   { cgi_mountains = Channel.fromPath([params.cgi_mountains_file, params.cgi_mountains_file + '.tbi'], checkIfExists: true).collect() } else { cgi_mountains = Channel.of([],[]) }
if (params.phastconselem_file)   { phastconselem = Channel.fromPath([params.phastconselem_file, params.phastconselem_file + '.tbi'], checkIfExists: true).collect() } else { phastconselem = Channel.of([],[]) }
if (params.encode_tfbs_file)     { encode_tfbs = Channel.fromPath([params.encode_tfbs_file, params.encode_tfbs_file + '.tbi'], checkIfExists: true).collect() } else { encode_tfbs = Channel.of([],[]) }

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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    sample_ch   = INPUT_CHECK.out.ch_sample

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

    // Prepare an input channel of vcf with sample names 
    vcf_ch = MPILEUP_SNV_CALL.out.vcf_ch
    ch_vcf = vcf_ch.join(name_ch)

    if (params.runSNVAnnotation){ 
        SNV_ANNOTATION(
            ch_vcf, 
            ref, 
            kgenome, 
            dbsnpsnv, 
            localcontrolwgs,
            localcontrolwes, 
            gnomadgenomes, 
            gnomadexomes, 
            annodb, 
            repeatmasker, 
            dacblacklist, 
            dukeexcluded, 
            hiseqdepth, 
            selfchain, 
            mapability, 
            simpletandemrepeats, 
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
            mirnas_sncrnas, 
            chr_prefix
        )
        ch_versions = ch_versions.mix(SNV_ANNOTATION.out.versions)

        //
        // SUBWORKFLOW: FILTER_SNVS: Filters SNVs
        //
        // input_ch= meta, annotated vcf, index, altbasequal, refbasequal, altreadpos, refreadpos, 
                    //sequence_spesific_error_plot_1, sequencing_spesific_error_plot_1, sequence_spesific_error_plot_2
                    //sequencing_spesific_error_plot_2, base_score_distribution_plot_1, base_score_distribution_plot_2
        vcf_ch   = SNV_ANNOTATION.out.vcf_ch
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
        }
        else{
            println "Skipping SNV filtering"
        }
    }
    else{
        println "Skipping SNV annotation and filtering"
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
