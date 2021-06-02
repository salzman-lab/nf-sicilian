//
// Check input samplesheet and get read channels
//
// Cribbed from https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/input_check.nf

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    // take:
    // samplesheet // file: /path/to/samplesheet.csv

    main:

    if (params.input_csv) { 
        ch_input_csv = file(params.input_csv) 
    } else if ( !(params.input || params.input_paths) ) { 
        exit 1, 'No input data specified with --input or --input_csv. Exiting!' 
    }

    star_output_params       = params.star_bam       && params.star_sj_out_tab       && params.reads_per_gene            && params.star_chimeric_junction
    star_output_params_paths = params.star_bam_paths && params.star_sj_out_tab_paths && params.star_reads_per_gene_paths && params.star_chimeric_junction_paths

    run_align = ! (star_output_params || star_output_params_paths)
    run_glm = !( params.sicilian_glm_output_paths || params.sicilian_glm_output )
    run_class_input = !( params.sicilian_class_input_paths || params.sicilian_class_input )

    /*
    * Initialize channels as empty
    */
    ch_reads          = Channel.empty()
    ch_input_csv      = Channel.empty()
    // STAR output
    ch_bam            = Channel.empty()
    ch_sj_out_tab     = Channel.empty()
    ch_reads_per_gene = Channel.empty()
    chimeric_junction = Channel.empty()
    // SICILIAN output
    ch_class_input    = Channel.empty()
    ch_glm_output     = Channel.empty()

    /*
    * Create a channel for input read files
    */
    single_end = params.single_end.toBoolean()
    if (run_align) {
        // No star
        if (params.input_paths) {
            if (params.single_end) {
                ch_reads = Channel.from(params.input_paths)
                    .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
                    .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
                    .map { it -> create_fastq_channels_from_filepairs( it, single_end, params.stranded) }
            } else {
                ch_reads = Channel.from(params.input_paths)
                    .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true), file(row[1][1], checkIfExists: true) ] ] }
                    .ifEmpty { exit 1, 'params.input_paths was empty - no input files supplied' }
                    // .map { it -> println "${it[0]}: it.getClass(): ${it.getClass()}" }
                    // .map { it -> println "${it[0]}: it[0].getClass(): ${it[0].getClass()}" }
                    // .map { it -> println "${it[0]}: it[1].getClass(): ${it[1].getClass()}" }
                    // .map { it -> println "${it[0]}: it[1][1].getClass(): ${it[1][1].getClass()}" }
                    .map { it -> create_fastq_channels_from_filepairs( it, single_end, params.stranded) }
            }
        } else if (params.input) {
            ch_reads = Channel.fromFilePairs(params.input, size: params.single_end ? 1 : 2)
                .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
                .map { it -> create_fastq_channels_from_filepairs( it, single_end, params.stranded) }
        } else if (params.input_csv) {
            samplesheet = file(params.input_csv, checkIfExists: true) 
            SAMPLESHEET_CHECK ( samplesheet )
                .splitCsv ( header:true, sep:',' )
                .map { create_fastq_channels(it) }
                .set { ch_reads }
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
                .dump ( tag: 'ch_class_input' )
        }
        if (params.sicilian_glm_output_paths) {
            ch_glm_output = Channel.from(params.sicilian_glm_output_paths)
                .map { row -> [ row[0], [ file(row[1][0], checkIfExists: true) ] ] }
                .ifEmpty { exit 1, 'params.sicilian_glm_output_paths was empty - no input files supplied' }
                .dump ( tag: 'ch_glm_output' )
        }
    }

    ch_reads.dump( tag: 'ch_reads' )
    ch_bam.dump ( tag: 'ch_bam' )
    ch_sj_out_tab.dump ( tag: 'ch_sj_out_tab' )
    ch_reads_per_gene.dump ( tag: 'ch_reads_per_gene' )
    ch_chimeric_junction.dump ( tag: 'ch_chimeric_junction' )
    ch_class_input.dump ( tag: 'ch_class_input' )
    ch_glm_output.dump ( tag: 'ch_glm_output' )

    run_align.dump( tag: 'run_align' )
    run_class_input.dump( tag: 'run_class_input' )
    run_glm.dump( tag: 'run_glm' )



    emit:
    reads                =  ch_reads        // channel: [ val(meta), [ reads ] ]
    // STAR output
    bam                  = ch_bam
    sj_out_tab           = ch_sj_out_tab
    reads_per_gene       = ch_reads_per_gene
    ch_chimeric_junction = chimeric_junction 
    // SICILIAN output
    class_input          = ch_class_input
    glm_output           = ch_glm_output
    // Booleans
    run_align            = run_align
    run_class_input      = run_class_input
    run_glm              = run_glm
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(LinkedHashMap row) {
    def meta = [:]
    meta.id               = row.sample_id
    meta.single_end       = row.single_end.toBoolean()
    meta.strandedness     = row.strandedness
    meta.concatenation_id = row.concatenation_id

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


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels_from_filepairs(ArrayList row, Boolean single_end, Boolean stranded) {
    def meta = [:]
    meta.id               = row[0]
    meta.single_end       = single_end
    meta.strandedness     = stranded
    meta.concatenation_id = meta.id

    // 1-indexed item is the reads
    def array = [ meta, row[1] ]
    return array
}