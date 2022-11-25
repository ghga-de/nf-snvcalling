/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSnvcalling.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

//
// Check mandatory parameters
//

if (params.input)         { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.annovar_path)  { annodb = Channel.fromPath(params.annovar_path + '/humandb/', checkIfExists: true ) } else { annodb = Channel.empty() }

// Set up reference depending on the genome choice
// NOTE: link will be defined by aoutomatic reference generation when the pipeline ready!
if (params.ref_type)
    {
    if (params.ref_type == 'hg37')
        { 
        def fa_file  = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/1KGRef_Phix/hs37d5_PhiX.fa"
        ref          = Channel.fromPath([fa_file,fa_file +'.fai'], checkIfExists: true).collect() 
        def chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hs37d5.fa.chrLenOnlyACGT_realChromosomes.tab'
        chrlength    = Channel.fromPath(chr_file, checkIfExists: true)
        chr_prefix   = Channel.value("")
        interval_ch  = Channel.fromList(['X','Y']) 
        }
    if (params.ref_type == 'hg19') 
        { 
        def fa_file  = "/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/sequence/hg19_chr/hg19_1-22_X_Y_M.fa"
        ref          = Channel.fromPath([fa_file,fa_file +'.fai'], checkIfExists: true).collect()
        def chr_file = '/omics/odcf/reference_data/legacy/ngs_share/assemblies/hg19_GRCh37_1000genomes/stats/hg19_1-22_X_Y_M.fa.chrLenOnlyACGT.tab'
        chrlength    = Channel.fromPath(chr_file, checkIfExists: true)
        chr_prefix   = Channel.value("chr")
        interval_ch  = Channel.value('chr1','chr2','chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY')  
        }
    }
else
{
    if (params.reference)      { ref = Channel.fromPath([params.reference,params.reference +'.fai'], checkIfExists: true).collect() } else { exit 1, 'Input reference file does not exist' }
    if (params.chrlength_file) { chrlength = Channel.fromPath(params.chrlength_file, checkIfExists: true) } else { exit 1, 'Chromosome length file does not exist'  }
    if (params.chr_prefix)     {chr_prefix= Channel.of(params.chr_prefix)} else {chr_prefix= Channel.of("")}
    interval_ch=Channel.of('1','2','3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y') 

}
// Annotation databases
if (params.k_genome)             { kgenome = Channel.fromPath([params.k_genome,params.k_genome +'.tbi'], checkIfExists: true).collect() } else { kgenome = Channel.empty() }
if (params.dbsnp_indel)          { dbsnpindel = Channel.fromPath([params.dbsnp_indel, params.dbsnp_indel + '.tbi'], checkIfExists: true).collect() } else { dbsnpindel = Channel.empty() }
if (params.dbsnp_snv)            { dbsnpsnv = Channel.fromPath([params.dbsnp_snv,params.dbsnp_snv +'.tbi' ], checkIfExists: true).collect() } else { dbsnpsnv = Channel.empty() }
if (params.exac_file)            { exac = Channel.fromPath([params.exac_file, params.exac_file + '.tbi'], checkIfExists: true).collect() } else { exac = Channel.empty() }
if (params.evs_file)             { evs = Channel.fromPath([params.evs_file, params.evs_file + '.tbi'], checkIfExists: true).collect() } else { evs = Channel.empty() }
if (params.local_control_wgs)    { localcontrolwgs = Channel.fromPath([params.local_control_wgs,params.local_control_wgs + '.tbi' ], checkIfExists: true).collect() } else { localcontrolwgs = Channel.empty() }
if (params.local_control_wes)    { localcontrolwes = Channel.fromPath([params.local_control_wes, params.local_control_wes + '.tbi'], checkIfExists: true).collect() } else { localcontrolwes = Channel.empty() }
if (params.gnomad_genomes)       { gnomadgenomes = Channel.fromPath([params.gnomad_genomes, params.gnomad_genomes + '.tbi'], checkIfExists: true).collect() } else { gnomadgenomes = Channel.empty() }
if (params.gnomad_exomes)        { gnomadexomes = Channel.fromPath([params.gnomad_exomes, params.gnomad_exomes + '.tbi'], checkIfExists: true).collect() } else { gnomadexomes = Channel.empty() }

// Annovar table folder
if (params.annovar_path)         { annodb = Channel.fromPath(params.annovar_path + '/humandb/', checkIfExists: true ) } else { annodb = Channel.empty() }
// Realiability files
if (params.repeat_masker)        { repeatmasker = Channel.fromPath([params.repeat_masker, params.repeat_masker + '.tbi'], checkIfExists: true).collect() } else { repeatmasker = Channel.empty() }
if (params.dac_blacklist)        { dacblacklist = Channel.fromPath([params.dac_blacklist, params.dac_blacklist + '.tbi'], checkIfExists: true).collect() } else { dacblacklist = Channel.empty() }
if (params.duke_excluded)        { dukeexcluded = Channel.fromPath([params.duke_excluded, params.duke_excluded + '.tbi'], checkIfExists: true).collect() } else { dukeexcluded = Channel.empty() }
if (params.hiseq_depth)          { hiseqdepth = Channel.fromPath([params.hiseq_depth, params.hiseq_depth + '.tbi'], checkIfExists: true).collect() } else { hiseqdepth = Channel.empty() }
if (params.self_chain)           { selfchain = Channel.fromPath([params.self_chain, params.self_chain + '.tbi'], checkIfExists: true).collect() } else { selfchain = Channel.empty() }
if (params.mapability_file)      { mapability = Channel.fromPath([params.mapability_file, params.mapability_file + '.tbi'], checkIfExists: true).collect() } else { mapability = Channel.empty() }
if (params.simple_tandemrepeats) { simpletandemrepeats = Channel.fromPath([params.simple_tandemrepeats, params.simple_tandemrepeats + '.tbi'], checkIfExists: true).collect() } else { simpletandemrepeats = Channel.empty() }
// Indel Deep Annotation files
if (params.enchancer_file)       { enchangers = Channel.fromPath([params.enchancer_file, params.enchancer_file + '.tbi'], checkIfExists: true).collect() } else { enchangers = Channel.of([],[]) }
if (params.cpgislands_file)      { cpgislands = Channel.fromPath([params.cpgislands_file, params.cpgislands_file + '.tbi'], checkIfExists: true).collect() } else { cpgislands = Channel.of([],[]) }
if (params.tfbscons_file)        { tfbscons = Channel.fromPath([params.tfbscons_file, params.tfbscons_file + '.tbi'], checkIfExists: true).collect() } else { tfbscons = Channel.of([],[]) }
if (params.encode_dnase_file)    { encode_dnase = Channel.fromPath([params.encode_dnase_file, params.encode_dnase_file + '.tbi'], checkIfExists: true).collect() } else { encode_dnase = Channel.of([],[]) }
if (params.mirnas_snornas_file)  { mirnas_snornas = Channel.fromPath([params.mirnas_snornas_file, params.mirnas_snornas_file + '.tbi'], checkIfExists: true).collect() } else { mirnas_snornas = Channel.of([],[]) }
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
    sample_ch = INPUT_CHECK.out.ch_sample

    //
    // MODULE: Extract sample name from BAM
    //
    GREP_SAMPLENAME(
        sample_ch
    )
    ch_versions = ch_versions.mix(GREP_SAMPLENAME.out.versions)

    // Prepare an input channel of sample with sample names
    name_ch=GREP_SAMPLENAME.out.samplenames
    ch_sample=sample_ch.join(name_ch)

    //
    // SUBWORKFLOW: MPILEUP_SNV_CALL: Call SNVs
    //
    
    MPILEUP_SNV_CALL(
        ch_sample, ref, interval_ch
    )
    ch_versions = ch_versions.mix(MPILEUP_SNV_CALL.out.versions)

    //
    // SUBWORKFLOW: SNV_ANNOTATION: Annotate SNVs
    //

    // Prepare an input channel of vcf with sample names 
    vcf_ch=MPILEUP_SNV_CALL.out.vcf_ch
    ch_vcf=vcf_ch.join(name_ch)

    if (params.runSNVAnnotation){ 
        SNV_ANNOTATION(
            ch_vcf,kgenome, dbsnpindel, exac, evs, localcontrolwgs,
        localcontrolwes, gnomadgenomes, gnomadexomes, annodb, repeatmasker, dacblacklist,
        dukeexcluded, hiseqdepth, selfchain, mapability, simpletandemrepeats, enchangers,
        cpgislands, tfbscons, encode_dnase, mirnas_snornas, cosmic, mirbase, mir_targets,
        cgi_mountains, phastconselem, encode_tfbs, chr_prefix
        )
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

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
