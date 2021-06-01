////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    params.input, 
    // params.multiqc_config,
    params.gtf, 
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


/*
========================================================================================
    SICILIAN-SPECIFIC FILES
========================================================================================
*/
//
// Create channel for domain file
//
ch_domain = file(params.domain, checkIfExists: true)


/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    SET MODULE PARAMETERS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()
def umitools_whitelist_options = modules['umitools_whitelist']
umitools_whitelist_options.args  += params.umitools_bc_pattern     ? " --bc-pattern='${params.umitools_bc_pattern}'"       : ''

def umitools_extract_options    = modules['umitools_extract']
umitools_extract_options.args  += params.umitools_bc_pattern     ? Utils.joinModuleArgs(["--bc-pattern='${params.umitools_bc_pattern}'"])       : ''
if (params.save_umi_intermeds)  { umitools_extract_options.publish_files.put('fastq.gz','') }

def star_genomegenerate_options = modules['star_genomegenerate']
if (!params.save_reference)     { star_genomegenerate_options['publish_files'] = false }

def gffread_options         = modules['gffread']
if (!params.save_reference) { gffread_options['publish_files'] = false }

def sicilian_createannotator_options = modules['sicilian_createannotator']

def star_align_options               = modules['star_align']

def sicilian_classinput_options      = modules['sicilian_classinput']
sicilian_classinput_options.args    += params.tenx ? Utils.joinModuleArgs(['--UMI_bar']) : ''

def sicilian_glm_options             = modules['sicilian_glm']
def sicilian_annsplices_options      = modules['sicilian_annsplices']
def sicilian_consolidate_options     = modules['sicilian_consolidate']
def sicilian_process_ci_10x_options  = modules['sicilian_process_ci_10x']
def sicilian_postprocess_options     = modules['sicilian_postprocess']

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

def multiqc_options         = modules['multiqc']
multiqc_options.args       += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { INPUT_CHECK              } from './subworkflows/local/input_check'                    addParams( options: [:] )
include { UMITOOLS_WHITELIST       } from './modules/local/umitools_whitelist'                  addParams( options: umitools_whitelist_options )
include { GET_SOFTWARE_VERSIONS    } from './modules/local/get_software_versions'               addParams( options: [publish_files : ['csv':'']]                      )
include { PREPARE_GENOME           } from './subworkflows/local/prepare_genome'                 addParams( 
    genome_options: publish_genome_options, 
    index_options: publish_index_options, 
    gffread_options: gffread_options,  
    star_index_options: star_genomegenerate_options, 
    sicilian_createannotator_options: sicilian_createannotator_options )
include { SICILIAN_CLASSINPUT } from './modules/local/sicilian/classinput.nf'        addParams( options: sicilian_classinput_options )
include { SICILIAN_GLM        } from './modules/local/sicilian/glm.nf'               addParams( options: sicilian_glm_options )
include { SICILIAN_ANNSPLICES } from './modules/local/sicilian/annsplices.nf'        addParams( options: sicilian_annsplices_options )

// Postprocessing of SICILIAN output, to consolidate output from mutliple samples into one summary file
include { SICILIAN_CONSOLIDATE    } from './modules/local/sicilian/consolidate.nf'  addParams( options: sicilian_consolidate_options )
include { SICILIAN_PROCESS_CI_10X } from './modules/local/sicilian/processci10x.nf' addParams( options: sicilian_process_ci_10x_options )
include { SICILIAN_POSTPROCESS    } from './modules/local/sicilian/postprocess.nf'  addParams( options: sicilian_postprocess_options )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC          } from './modules/nf-core/software/multiqc/main.nf'                     addParams( options: multiqc_options                                   )
include { UMITOOLS_EXTRACT } from './modules/nf-core/software/umitools/extract/main.nf' addParams( options: umitools_extract_options )
include { STAR_ALIGN       } from './modules/nf-core/software/star/align/main.nf' addParams( options: star_align_options )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report      = []
def pass_percent_mapped = [:]
def fail_percent_mapped = [:]

workflow SICILIAN {

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    // PREPARE_GENOME (
    //     prepareToolIndices
    // )
    ch_software_versions = Channel.empty()
    // ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK ()
    ch_reads        = INPUT_CHECK.out.reads
    run_align       = INPUT_CHECK.out.run_align
    run_class_input = INPUT_CHECK.out.run_class_input
    run_glm         = INPUT_CHECK.out.run_glm
    ch_reads.dump( tag: 'ch_reads' )


    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME ()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.star_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))

    ch_star_multiqc  = Channel.empty()
    if (run_align) {
        if (!params.skip_umitools) {
        /*
        * MODULE: Create a whitelist of UMIs from the data
        */
        UMITOOLS_WHITELIST ( 
            ch_reads,
        )
        ch_software_versions = ch_software_versions.mix(UMITOOLS_WHITELIST.out.version.ifEmpty(null))


        UMITOOLS_EXTRACT ( ch_reads, UMITOOLS_WHITELIST.out.whitelist ).reads.set { umi_reads }
        ch_software_versions = ch_software_versions.mix(UMITOOLS_EXTRACT.out.version.ifEmpty(null))
        } else {
            umi_reads = ch_reads
        }

        STAR_ALIGN (
            umi_reads,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
        )
        ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.first().ifEmpty(null))
        ch_bam               = STAR_ALIGN.out.bam
        ch_sj_out_tab        = STAR_ALIGN.out.sj_out_tab
        ch_chimeric_junction = STAR_ALIGN.out.chimeric_out_junction
        ch_reads_per_gene    = STAR_ALIGN.out.reads_per_gene
        ch_star_multiqc      = STAR_ALIGN.out.log_final
    } else {
        // Skipping alignment, because already have all the files from the input
        ch_bam = INPUT_CHECK.out.bam
        ch_sj_out_tab = INPUT_CHECK.out.sj_out_tab
        ch_chimeric_junction = INPUT_CHECK.out.chimeric_out_junction
        ch_reads_per_gene = INPUT_CHECK.out.reads_per_gene
    }



    if (run_class_input) {
        SICILIAN_CLASSINPUT (
            ch_bam,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.sicilian_annotator,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_CLASSINPUT.out.version.ifEmpty(null))
        ch_class_input = SICILIAN_CLASSINPUT.out.class_input
    } else {
        ch_class_input = INPUT_CHECK.out.class_input
    }

    if (run_glm) {
        ch_classinput_sjouttab_chimericjunctions_readspergene = ch_class_input.join(
            ch_sj_out_tab
        ).join( 
            ch_chimeric_junction
        ).join(
            ch_reads_per_gene
        )

        SICILIAN_GLM (
            PREPARE_GENOME.out.gtf,
            ch_domain,
            PREPARE_GENOME.out.sicilian_exon_bounds,
            PREPARE_GENOME.out.sicilian_splices,
            ch_classinput_sjouttab_chimericjunctions_readspergene
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_GLM.out.version.ifEmpty(null))
        ch_glm_output = SICILIAN_GLM.out.glm_output

        SICILIAN_ANNSPLICES(
            SICILIAN_GLM.out.sicilian_called_splices,
            PREPARE_GENOME.out.sicilian_exon_bounds,
            PREPARE_GENOME.out.sicilian_splices,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_GLM.out.version.ifEmpty(null))
    } else {
        ch_glm_output = INPUT_CHECK.out.glm_output
    }

    SICILIAN_CONSOLIDATE(
        // Take the 2nd (1-index) item, which is the file only, and not the sample id
        ch_glm_output.collect{ it[1] }.dump(tag: 'ch_glm_output_collected')
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_CONSOLIDATE.out.version.ifEmpty(null))

    ch_class_input
        .dump( tag: 'class_input' )
        .multiMap { 
            meta, class_input ->
                concatenation_ids: meta.concatenation_id
                class_inputs: class_input
        }
        .set { class_input_multimap }
    ch_process_ci_concatenation_ids = class_input_multimap.concatenation_ids
        .collect()
        .dump( tag: 'class_input_concatenation_ids_collected' )
    ch_process_ci_class_inputs = class_input_multimap.class_inputs
        .collect()
        .dump( tag: 'class_input_files_collected' )

    SICILIAN_PROCESS_CI_10X (
        ch_process_ci_concatenation_ids,
        ch_process_ci_class_inputs,
        PREPARE_GENOME.out.gtf,
        PREPARE_GENOME.out.sicilian_exon_bounds,
        PREPARE_GENOME.out.sicilian_splices,
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_PROCESS_CI_10X.out.version.ifEmpty(null))


    SICILIAN_POSTPROCESS (
        SICILIAN_PROCESS_CI_10X.out.sicilian_junctions_tsv,
        SICILIAN_CONSOLIDATE.out.glm_consolidated
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_POSTPROCESS.out.version.ifEmpty(null))


    ch_software_versions
        .flatten()
        .dump(tag: 'ch_software_versions__flatten')
        .map { it -> if (it) [ it.baseName, it ] }
        .dump(tag: 'ch_software_versions__flatten__map')
        .groupTuple()
        .dump(tag: 'ch_software_versions__flatten__map__grouptuple')
        .map { it[1][0] }
        .dump(tag: 'ch_software_versions__flatten__map__grouptuple__map')
        .collect()
        .dump(tag: 'ch_software_versions__flatten__map__grouptuple__map__collect')
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions
    )

    /*
     * MultiQC
     */
    if (!params.skip_multiqc) {
        workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            ch_star_multiqc.collect{it[1]}.ifEmpty([]),
        )
        multiqc_report = MULTIQC.out.report.toList()
    }
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report, fail_percent_mapped)
    Completion.summary(workflow, params, log, fail_percent_mapped, pass_percent_mapped)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(Tuple row, Boolean single_end) {
    def meta = [:]
    meta.id           = row[0]
    meta.single_end   = single_end.toBoolean()
    meta.strandedness = row.strandedness

    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}