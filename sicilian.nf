////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Check input path parameters to see if they exist
checkPathParamList = [
    // params.input, 
    // params.multiqc_config,
    params.gtf, 
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

star_output_params = params.star_sj_out_tab && params.reads_per_gene && params.star_chimeric_junction && params.star_bam
star_output_params_paths = params.star_bam_paths && params.star_sj_out_tab_paths && params.star_reads_per_gene_paths && params.star_chimeric_junction_paths

run_align = ! (star_output_params || star_output_params_paths)
run_glm = !( params.sicilian_glm_output_paths || params.sicilian_glm_output )
run_class_input = !( params.sicilian_class_input_paths || params.sicilian_class_input )

/*
 * Create a channel for input read files
 */
if (run_align) {
    // No star
    if (params.input_paths) {
        if (params.single_end) {
            ch_reads = Channel.from(params.input_paths)
                .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
                .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
        } else {
            ch_reads = Channel.from(params.input_paths)
                .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
                .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
        }
    } else {
        ch_reads = Channel.fromFilePairs(params.input, size: params.single_end ? 1 : 2)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
    }
} else {
    if (params.star_bam_paths) {
        ch_bam = Channel.from(params.star_bam_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.star_bam_paths was empty - no input files supplied' }
    }   
    if (params.star_sj_out_tab_paths) {
        ch_sj_out_tab = Channel.from(params.star_sj_out_tab_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.star_sj_out_tab_paths was empty - no input files supplied' }
    }
    if (params.star_reads_per_gene_paths) {
        ch_reads_per_gene = Channel.from(params.star_reads_per_gene_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.star_reads_per_gene_paths was empty - no input files supplied' }
    }
    if (params.star_chimeric_junction_paths) {
        ch_chimeric_junction = Channel.from(params.star_chimeric_junction_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.star_chimeric_junction_paths was empty - no input files supplied' }
    }
    if (params.sicilian_class_input_paths) {
        ch_class_input = Channel.from(params.sicilian_class_input_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.ch_class_input was empty - no input files supplied' }
    }
    if (params.sicilian_glm_output_paths) {
        ch_glm_output = Channel.from(params.sicilian_glm_output_paths)
            .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
            .ifEmpty { exit 1, 'params.sicilian_glm_output_paths was empty - no input files supplied' }
    }
}


//
// Create channel for domain file
//
ch_domain = file(params.domain, checkIfExists: true)


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

def star_align_options            = modules['star_align']
def sicilian_createannotator_options            = modules['sicilian_createannotator']

def sicilian_classinput_options    = modules['sicilian_classinput']
sicilian_classinput_options.args   += params.tenx ? Utils.joinModuleArgs(['--UMI_bar']) : ''

def sicilian_glm_options    = modules['sicilian_glm']
def sicilian_annsplices_options    = modules['sicilian_annsplices']
def sicilian_consolidate_options    = modules['sicilian_consolidate']
def sicilian_process_ci_10x_options    = modules['sicilian_process_ci_10x']
def sicilian_postprocess_options    = modules['sicilian_postprocess']

def publish_genome_options = params.save_reference ? [publish_dir: 'genome']       : [publish_files: false]
def publish_index_options  = params.save_reference ? [publish_dir: 'genome/index'] : [publish_files: false]

// Import modules
include { UMITOOLS_WHITELIST       } from './modules/local/umitools_whitelist'          addParams( options: umitools_whitelist_options )
include { UMITOOLS_EXTRACT         } from './modules/nf-core/software/umitools/extract/main.nf'   addParams( options: umitools_extract_options )
include { GET_SOFTWARE_VERSIONS    } from './modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']]                      )
include { PREPARE_GENOME           } from './subworkflows/local/prepare_genome' addParams( 
    genome_options: publish_genome_options, 
    index_options: publish_index_options, 
    gffread_options: gffread_options,  
    star_index_options: star_genomegenerate_options, 
    sicilian_createannotator_options: sicilian_createannotator_options )
include { STAR_ALIGN          } from './modules/nf-core/software/star/align/main.nf'          addParams( options: star_align_options )
include { SICILIAN_CLASSINPUT } from './modules/local/sicilian/classinput.nf'          addParams( options: sicilian_classinput_options )
include { SICILIAN_GLM        } from './modules/local/sicilian/glm.nf'          addParams( options: sicilian_glm_options )
include { SICILIAN_ANNSPLICES } from './modules/local/sicilian/annsplices.nf'          addParams( options: sicilian_annsplices_options )

// Postprocessing of SICILIAN output, to consolidate output from mutliple samples into one summary file
include { SICILIAN_CONSOLIDATE    } from './modules/local/sicilian/consolidate.nf'          addParams( options: sicilian_consolidate_options )
include { SICILIAN_PROCESS_CI_10X } from './modules/local/sicilian/processci10x.nf'          addParams( options: sicilian_process_ci_10x_options )
include { SICILIAN_POSTPROCESS    } from './modules/local/sicilian/postprocess.nf'          addParams( options: sicilian_postprocess_options )

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

    // TODO: Add INPUT_CHECK subworkflow to allow for samplesheet input
    // Example usage: https://github.com/nf-core/rnaseq/blob/0fcbb0ac491ecb8a80ef879c4f3dad5f869021f9/workflows/rnaseq.nf#L250
    // Subwokflow: https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/input_check.nf

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    PREPARE_GENOME ()
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.star_version.ifEmpty(null))
    ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.gffread_version.ifEmpty(null))


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
        ch_bam = STAR_ALIGN.out.bam
        ch_sj_out_tab = STAR_ALIGN.out.sj_out_tab
        ch_chimeric_junction = STAR_ALIGN.out.chimeric_out_junction
        ch_reads_per_gene = STAR_ALIGN.out.reads_per_gene
    }



    if (run_class_input) {
        SICILIAN_CLASSINPUT (
            ch_bam,
            PREPARE_GENOME.out.gtf,
            PREPARE_GENOME.out.sicilian_annotator,
        )
        ch_software_versions = ch_software_versions.mix(SICILIAN_CLASSINPUT.out.version.ifEmpty(null))
        ch_class_input = SICILIAN_CLASSINPUT.out.class_input
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
    }

    SICILIAN_CONSOLIDATE(
        // Take the 2nd (1-index) item, which is the file only, and not the sample id
        ch_glm_output.collect{ it[1] }.dump(tag: 'ch_glm_output_collected')
    )
    ch_software_versions = ch_software_versions.mix(SICILIAN_CONSOLIDATE.out.version.ifEmpty(null))

    ch_class_input
        .dump( tag: 'class_input' )
        .multiMap {
            sample_ids: it[0]
            class_inputs: it[1]
        }
        .set { class_input_multimap }
    ch_process_ci_sample_ids = class_input_multimap.sample_ids.collect().dump( tag: 'class_input_sample_ids_collected')
    ch_process_ci_class_inputs = class_input_multimap.class_inputs.collect().dump( tag: 'class_input_files_collected')

    SICILIAN_PROCESS_CI_10X (
        ch_process_ci_sample_ids,
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
/// STAR needs  --sjdbGTFfile {} .format(gtf_file) option


    // /*
    //  * SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //  */
    // ch_genome_bam                 = Channel.empty()
    // ch_genome_bai                 = Channel.empty()
    // ch_samtools_stats             = Channel.empty()
    // ch_samtools_flagstat          = Channel.empty()
    // ch_samtools_idxstats          = Channel.empty()
    // ch_star_multiqc               = Channel.empty()
    // ch_aligner_pca_multiqc        = Channel.empty()
    // ch_aligner_clustering_multiqc = Channel.empty()
    // ALIGN_STAR (
    //     ch_trimmed_reads,
    //     PREPARE_GENOME.out.star_index,
    //     PREPARE_GENOME.out.gtf
    // )
    // ch_genome_bam        = ALIGN_STAR.out.bam
    // ch_genome_bai        = ALIGN_STAR.out.bai
    // ch_samtools_stats    = ALIGN_STAR.out.stats
    // ch_samtools_flagstat = ALIGN_STAR.out.flagstat
    // ch_samtools_idxstats = ALIGN_STAR.out.idxstats
    // ch_star_multiqc      = ALIGN_STAR.out.log_final
    // ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.star_version.first().ifEmpty(null))
    // ch_software_versions = ch_software_versions.mix(ALIGN_STAR.out.samtools_version.first().ifEmpty(null))

    // /*
    //     * SUBWORKFLOW: Count reads from BAM alignments using Salmon
    //     */
    // QUANTIFY_STAR_SALMON (
    //     ALIGN_STAR.out.bam_transcript,
    //     ch_dummy_file,
    //     PREPARE_GENOME.out.transcript_fasta,
    //     PREPARE_GENOME.out.gtf,
    //     true
    // )
    // ch_software_versions = ch_software_versions.mix(QUANTIFY_STAR_SALMON.out.salmon_version.first().ifEmpty(null))
    // ch_software_versions = ch_software_versions.mix(QUANTIFY_STAR_SALMON.out.tximeta_version.first().ifEmpty(null))
    // ch_software_versions = ch_software_versions.mix(QUANTIFY_STAR_SALMON.out.summarizedexperiment_version.ifEmpty(null))

    // if (!params.skip_qc & !params.skip_deseq2_qc) {
    //     DESEQ2_QC_STAR_SALMON (
    //         QUANTIFY_STAR_SALMON.out.merged_counts_gene_length_scaled,
    //         ch_pca_header_multiqc,
    //         ch_clustering_header_multiqc
    //     )
    //     ch_aligner_pca_multiqc        = DESEQ2_QC_STAR_SALMON.out.pca_multiqc
    //     ch_aligner_clustering_multiqc = DESEQ2_QC_STAR_SALMON.out.dists_multiqc
    //     ch_software_versions          = ch_software_versions.mix(DESEQ2_QC_STAR_SALMON.out.version.ifEmpty(null))
    // }


    /*
     * MODULE: Pipeline reporting
     */
    // GET_SOFTWARE_VERSIONS ( 
    //     ch_software_versions.map { it }.collect()
    // )

    /*
     * MultiQC
     */
    // if (!params.skip_multiqc) {
    //     workflow_summary    = Schema.params_summary_multiqc(workflow, params.summary_params)
    //     ch_workflow_summary = Channel.value(workflow_summary)

    //     MULTIQC (
    //         ch_multiqc_config,
    //         ch_multiqc_custom_config.collect().ifEmpty([]),
    //         GET_SOFTWARE_VERSIONS.out.yaml.collect(),
    //         ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
    //         // ch_fail_mapping_multiqc.ifEmpty([]),
    //         // ch_fail_strand_multiqc.ifEmpty([]),
    //         // FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
    //         // FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
    //         // FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
    //         // ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_star_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_hisat2_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_rsem_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_salmon_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_samtools_stats.collect{it[1]}.ifEmpty([]),
    //         // ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
    //         // ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
    //         // ch_markduplicates_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_featurecounts_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_aligner_pca_multiqc.collect().ifEmpty([]),
    //         // ch_aligner_clustering_multiqc.collect().ifEmpty([]),
    //         // ch_pseudoaligner_pca_multiqc.collect().ifEmpty([]),
    //         // ch_pseudoaligner_clustering_multiqc.collect().ifEmpty([]),
    //         // ch_preseq_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_qualimap_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_dupradar_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_bamstat_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_inferexperiment_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_innerdistance_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_junctionannotation_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_junctionsaturation_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_readdistribution_multiqc.collect{it[1]}.ifEmpty([]),
    //         // ch_readduplication_multiqc.collect{it[1]}.ifEmpty([])
    //     )
    //     multiqc_report = MULTIQC.out.report.toList()
    // }
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